import numpy as np
import matplotlib.pyplot as plt

from scipy.spatial import cKDTree
from geom_primitives import hexagonal_lattice, periodic_domain
from geom_primitives import flat2torus, cartesian_lattice, circular_porous_domain

from geom_primitives import save_to_hdf5, distances, save2txt, \
    save_support

def normalize(v):
    """Normalize a vector.

    Args:
        v (array-like): A one-dimensional vector type object.

    Returns:
        The normalized vector.
    """
    norm=np.linalg.norm(v, ord=1)
    if norm==0:
        norm=np.finfo(v.dtype).eps
    return v/norm


def tiled_point_cloud(bbox,xy,field=None):
    """Construct a a periodic point cloud tiling.

    This function returns a larger point cloud, built as a tiling
    of a smaller one. The larger point cloud can then be passed to 
    a spatial search structure for periodic searching. The indexes
    in the original point cloud can be recovered by taking the
    remainder of the division between point cloud index and the number
    of points in the point cloud. Also see `find_periodic_supports`.

    TODO: Support for point clouds periodic in only one dimension.

    Args:
        bbox (list): Bounding box of the domain.
        xy (ndarray): Numpy array of x- and y-coordinates of an initial point cloud.
        field (ndarray, optional): A scalar field that we would like to tile.
            This is mainly used for visualization, say if we would like to plot
            our periodic field of points.

    Returns:
        If field is present returns a tuple (txy,tfield) if the tiled point cloud
        and corresponding field. If field is None, only the tiled point cloud is returned.
    """

    w = bbox[1] - bbox[0]
    h = bbox[3] - bbox[2]

    n = xy.shape[0]

    tiled_xy = np.empty((9*n,2))

    if field is not None:
        tiled_field = np.empty(9*n)

    # The default arrangement of tiles (D2Q9 numbering):
    #
    # x------x------x------x
    # | 6    | 2    | 5    |
    # |      |      |      |
    # x------x------x------x
    # | 3    | 0    | 1    |
    # |      |      |      |
    # x------x------x------x
    # | 7    | 4    | 8    |
    # |      |      |      |
    # x------x------x------x
    #
    # Tile 0 represents the original point cloud.
    # The other tiles are created by translating the original
    # coordinates according to the tile vector, and the 
    # bounding box data.

    tiles = [(0,0),(1,0),(0,1),(-1,0),(0,-1),(1,1),(-1,1),(-1,-1),(1,-1)]

    for i,tile in enumerate(tiles):

        tiled_xy[i*n:i*n+n,0] = xy[:,0] + w*tile[0] 
        tiled_xy[i*n:i*n+n,1] = xy[:,1] + h*tile[1]

        if field is not None:
            tiled_field[i*n:i*n+n] = field

    if field is not None:
        return tiled_xy, tiled_field
    else:
        return tiled_xy


def find_periodic_supports(bbox,xy,ss,centers=None,**kwargs):
    """Find support stencils for a fully periodic domain.

    This function constructs a tiled point cloud, which is used to
    populate a kD-tree for fast spatial searches. We use the
    remainder trick to find the coordinates in the original array.

    Args:
        bbox (list): Bounding box of the domain.
        xy (array): Array of two-dimensional Cartesian point coordinates.
        ss (int): Stencil size.
            The number of points included in each stencil.
        centers (array): Optional array of center locations

    Returns:
        Adjacency list of the periodic supports. 
    """

    n = xy.shape[0]
    temp_xy = tiled_point_cloud(bbox,xy)
    tree = cKDTree(temp_xy)

    if centers is not None:
        _dists, supports = tree.query(centers,ss)
    else:
        # Determine the supports of the points xy
        _dists, supports = tree.query(xy,ss)

    idxs = np.remainder(supports,n)

    ret_tree = kwargs.get('ret_tree',False)
    if (ret_tree):
        return idxs, tree
    else:
        return idxs

def dist(p1,p2,w,h):
    """Distance between two points in a doubly periodic AABB.

    Args:
        p1 (array-like): 2-dimensional vector
        p2 (array-like): 2-dimensional vector
        w (float): domain width
        h (float): domain height

    Returns:
        Distance between points p1 and p2.
    """

    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    if (dx > 0.5*w):
        dx -= w 
    if (dx <= -0.5*w):
        dx += w
    if (dy > 0.5*h):
        dy -= h
    if (dy <= -0.5*h):
        dy += h

    return np.hypot(dx,dy)

def points_to_bbox(bbox,pos):
    """Corrects a set of points to be within the periodic bounding box

    Moves the point coordinates back within the periodic bound box.

    Args:
        bbox (list): The periodic bounding box.
        pos (array-like): N times 2 array of point coordinates.

    Returns:
        The values in `pos` are modified in-place.
    """

    assert pos.shape[1] == 2

    delta = [bbox[1] - bbox[0], bbox[3] - bbox[2]]
    for i,b in enumerate([0,2]):

        m_sml = pos[:,i] < bbox[b]
        m_lrg  = pos[:,i] >= bbox[b+1]

        pos[m_sml,i] += delta[i]
        pos[m_lrg,i] -= delta[i]


def relax_nodes(bbox,xy,iters,which=None,ss=7,rfunc=None,heat=(1,1),po=2,rebuild_after=1):
    """Relax point-cloud according to a distance function.

    Args:
        bbox (list): Boundinx box of the periodic domain.
        xy (array-like): N times 2 array of coordinates.
        iters (int): Number of iterations
        which (list, optional): A subset of nodes to apply relaxation. 
            Default is all nodes. This can be useful if we only want to
            relax internal nodes, and keep the boundary nodes fixed.
        ss (int): Support size for the repelling force calculation.
        rfunc (function, optional): Distance function used during relaxation.
        heat (tuple, optional): Starting and final heat values for annealing.
            The default values are the same, meaning no annealing is performed.
        po (int, optional): Order of the force potential.
        rebuild_after (int, optional): Number of iterations upon which we rebuild the tree.
            Due to point movement, the initial spatial search tree will give
            wrong results. This is used to reinitialize the tree every `rebuild_after` steps.

    Returns:
        No return value. The points in `xy` are modified in-place.
    """

    w = bbox[1] - bbox[0]
    h = bbox[3] - bbox[2]

    n = xy.shape[0]

    if which is None:
        which = list(range(xy.shape[0]))

    for i in range(iters):
        
        temp_xy = tiled_point_cloud(bbox,xy)

        # Build KD tree of periodic point set
        if i == 0 or i % rebuild_after == 0:
            print("Rebuilding KD-Tree")
            tree = cKDTree(temp_xy)

        # Find supports
        _dists, supports = tree.query(temp_xy,ss)

        # supports = np.remainder(supports,n)

        # Linear annealing factor
        iheat, fheat = heat
        af = fheat + (iheat - fheat)*(iters - i)/iters

        # Initalize node displacement vector

        dpos = np.zeros_like(xy)

        for k, supp in zip(which,supports[which,:]):

            p = xy[k,:]
            # d0 = dist(p,xy[supp[1] % n,:],w,h)
            d0 = np.linalg.norm(p - temp_xy[supp[1],:])

            temp_supp = np.remainder(supp[1:],n)
            # print(temp_supp)
            # r = xy[temp_supp,:] - p
            r = temp_xy[supp[1:],:] - p
            rx = r[:,0]
            ry = r[:,1]
            rnorm = np.sqrt(rx**2 + ry**2)

            assert np.all(rnorm > 1e-15)

            if callable(rfunc):
                r0 = rfunc(xy[temp_supp,:])
            else:
                r0 = rfunc*np.ones_like(temp_supp,dtype=np.double)

            rnorm = (np.abs(r0)/rnorm)**po
            rnorm = rnorm[np.newaxis,:]
            dpos[k,:] = rnorm.dot(r)

            if not callable(rfunc):
                if rfunc < 0:
                    dnorm = normalize(dpos[k,:])
                    dpos[k,:] = 0.2*dnorm*d0*af

            dpos[k,:] *= -af*d0
        # Move all points
        xy += dpos
        points_to_bbox(bbox,xy)


def main(domain_size=99,support_size=21,solid_fraction=0.3,h5file="relaxed.h5",h5group="relaxed",
    bb_correction=False,display=False,**kwargs):
    """Construct a point cloud for the square array of cylinders benchmark
    """

    # (xy, _p, _b), dx = hexagonal_lattice(bbox,dx,retstep=True)
    # xy, _p, _b = periodic_domain(bbox,dx)

    if "seed" in kwargs:
        seed = kwargs["seed"]
    else:
        seed = None
    print("Seed = ", seed)

    bbox = [0,1,0,1]
    dx = 1./domain_size

    # Fraction of solid in porous cylinder
    radius = np.sqrt(solid_fraction/np.pi)
    
    # Make the radius bigger for half a step due to bounceback
    if bb_correction:
        radius += 0.5*dx

    print("solid_fraction = ", solid_fraction )

    center = (0.5,0.5)
    xy, t, bmap, normals = circular_porous_domain(center,radius,dx,seed=seed)


    xy *= (1./dx)
    bbox = (1./dx)*np.asarray(bbox)
    radius *= (1./dx)
    
    print("radius = ", radius)
    print("bbox = ", bbox)


    which = np.flatnonzero(bmap < 0)
    print(which)

    xy_old = np.copy(xy)


    # def rfunc(pos):
        # return dx


    if 'num_iters' in kwargs:
        num_iters = kwargs['num_iters']
    else:
        num_iters = 50

    if 'po' in kwargs:
        po = kwargs['po']
    else:
        po = 3

    rfunc = -1

    relax_nodes(bbox,xy,num_iters,which=which,ss=5,heat=(1,0),rfunc=rfunc,po=po,rebuild_after=1)
    
    ss = support_size
    supp = find_periodic_supports(bbox,xy,ss)

    save_to_hdf5(h5file,h5group,xy,t,bmap,
        normals=normals,
        bbox=bbox,
        supports=supp)

    if display:
        temp_xy = tiled_point_cloud(bbox,xy)

        # plt.scatter(xy_old[:,0],xy_old[:,1],marker='o')
        plt.scatter(temp_xy[:,0],temp_xy[:,1],marker='o')
        plt.scatter(xy[supp[121,:],0],xy[supp[121,:],1],marker='x')

        import matplotlib.patches as patches

        corner = [bbox[0],bbox[2]]
        w = bbox[1] - bbox[0]
        h = bbox[3] - bbox[2]
        
        tiles = [(0,0),(1,0),(0,1),(-1,0),(0,-1),(1,1),(-1,1),(-1,-1),(1,-1)]

        ax = plt.gca()
        for t in tiles:

            tcorner = np.asarray(corner) + np.asarray(t)*np.asarray([w,h])
            p = patches.Rectangle(tcorner,w,h,fill=False) 


            ax.add_patch(p)

        plt.axis('equal')

        plt.show()

def taylor_green_scattered(sz,ss,display=False,ret=False,**kwargs):

    bbox = [0.,1.,0.,1.]
    dx = 1./sz

    seed = kwargs.get('seed',2019)

    p,t,b = periodic_domain(bbox=bbox,dx=dx,seed=seed)

    iters = kwargs.get('iters',40)
    heat = kwargs.get('heat',(1.,1.))
    po = kwargs.get('po',4)

    rfunc = -1
    which = np.flatnonzero(b < 0)
    relax_nodes(bbox,p,iters,ss=7,rfunc=rfunc,heat=heat,po=po,rebuild_after=1)

    p *= (1./dx)
    bbox = (1./dx)*np.asarray(bbox)

    for s in ss:
        supp = find_periodic_supports(bbox,p,s)    

        normals = np.empty((0,2))

        groupname = "tg_scattered_{:d}_{:d}".format(sz,s)
        filename = groupname + ".h5"
        save_to_hdf5(filename,groupname,p,t,b,normals=normals,bbox=bbox,supports=supp)
        save2txt(groupname,p,supp)

    interp = kwargs.get('interp',False)
    if (interp):
        # Calculate interpolation points
        (cart_p,cart_t,cart_b) = cartesian_lattice(bbox=bbox,dx=(1.,1.),add_half=True)
        _dists, idxs = tree.query(cart_p,k=1)
        idxs = np.remainder(idxs,p.shape[0])
        with open(groupname + '_interpxy.txt','w') as ff:
            ff.write('{:d}\n'.format(idxs.shape[0]))
            for i, x, y in zip(idxs,cart_p[:,0], cart_p[:,1]):
                row = '{:d} {:27.14e} {:27.14e}\n'.format(i,x,y)
                ff.write(row)

    semilagrangian = kwargs.get('semilagrangian',False)
    if (semilagrangian):
        cx = [0,1,0,-1,0,1,-1,-1,1]
        cy = [0,0,1,0,-1,1,1,-1,-1]
        dt = kwargs.get('dt',1.0)

        supports = []
        for cxi, cyi in zip(cx,cy):

            # The lattice vector
            c = dt*np.array([cxi,cyi])

            # Subtract direction from center (automatic broadcasting)
            pshifted = p - c 

            # Find supports of the shifted centers
            supp = find_periodic_supports(bbox,p,s,pshifted)
            supports.append(supp)

        # Now we have 9 support graphs
        for i in range(9):
            save_support('{}_sl{:d}'.format(groupname,i),supports[i])


    if display:
        # from torus_playground import plot_torus
        # fig3d, ax3d = plot_torus(x,y,z)

        grid_lines = np.linspace(0,sz,sz+1)
        plt.scatter(p[:,0],p[:,1])
        plt.scatter(p[supp[400,:],0],p[supp[400,:],1],marker='x')
        plt.hlines(grid_lines,0,sz)
        plt.vlines(grid_lines,0,sz)

        plt.axis("equal")
        plt.show()

    if ret:
        return p, t, b, normals, supp, bbox


def taylor_green_cartesian(sz,ss,display=False,ret=False):

    bbox = [0.,1.,0.,1.]
    dx = 1./sz

    (p,t,b), (dx, dy) = cartesian_lattice(bbox=bbox,dx=(dx,dx),add_half=True,retstep=True)
    print("dx, dy = ", dx,dy)

    p *= (1./dx)
    bbox = (1./dx)*np.asarray(bbox)

    supp = find_periodic_supports(bbox,p,ss)

    normals = np.empty((0,2))

    groupname = "tg_cartesian_{}_{}".format(str(sz),str(ss))
    filename = groupname + ".h5"
    save_to_hdf5(filename,groupname,p,t,b,normals=normals,bbox=bbox,supports=supp)
    save2txt(groupname,p,supp)

    if display:
        # from torus_playground import plot_torus
        # fig3d, ax3d = plot_torus(x,y,z)

        grid_lines = np.linspace(0,sz,sz+1)
        plt.scatter(p[:,0],p[:,1])
        plt.scatter(p[supp[400,:],0],p[supp[400,:],1],marker='x')
        plt.hlines(grid_lines,0,sz)
        plt.vlines(grid_lines,0,sz)
        plt.axis("equal")
        plt.show()

    if ret:
        return p, t, b, normals, supp, bbox

def make_cartesian_grids(display=False):

    for sz in [32,64,128,256]:
        for ss in [9,13,21,25,29,37,45,49]:
            taylor_green_cartesian(sz, ss,display=display)


def make_scattered_grids(display=False,**kwargs):

    ss = [9,13,21,37]
    for sz in [64]:
        taylor_green_scattered(sz, ss,display=display,**kwargs)

def generate_grids_for_paper(test):

    seed = 2108 # 50 iters
    seed = 1009 # 100 iters
    seed = 1702 # 150 iters
    num_iters = 150
    po = 3
    support_size = 21

    if test == 1:
        name = "cylinder1"
        sfracs = [0.2,0.3,0.4,0.5]
        sizes = [33,66,132,264]

        for sz in sizes:
            for sf in sfracs:

                group = name + "_" + str(sz) + "_" + str(sf)
                fname = group + ".h5"

                main(domain_size=sz,
                     support_size=support_size,
                     solid_fraction=sf,
                     h5file=fname,
                     h5group=group,
                     num_iters=num_iters,
                     po=po,
                     seed=seed)
    elif test == 2:

        name = "cylinder_2"
        sfrac = 0.4
        sizes = [33 + i for i in range(0,103,3)]

        for sz in sizes:

            group = name + "_" + str(sz) + "_" + str(sfrac)
            fname = group + ".h5"

            main(domain_size=sz,
                 support_size=support_size,
                 solid_fraction=sfrac,
                 h5file=fname,
                 h5group=group,
                 num_iters=num_iters,
                 po=po,
                 seed=seed)



if __name__ == '__main__':


    # make_cartesian_grids(display=False)
    # make_scattered_grids(display=True,iters=200,heat=(1,0.5))
    
    # CFL = 1.6
    dt = 1.13137084989848

    # CFL = 6
    dt = 8.48528137423857

    make_scattered_grids(display=True,iters=200,heat=(1,0.5),
        semilagrangian=True,
        dt=dt)



    # taylor_green_scattered(32, 21,display=True,iters=100,heat=(1.0,0.2))
    # taylor_green_scattered(64, 21,display=True,iters=100,heat=(1.0,0.2))


    # rbf_test()

    # for t in [1,2]:
        # generate_grids_for_paper(t)
    
    # generate_grids_for_paper(2)

    # bb_correction = False
    # num_iters=100
    # seed = 2019
    # po=3

    # main(domain_size=1*33,support_size=35,solid_fraction=0.4,h5file="cylinder33.h5",h5group="cylinder33",num_iters=num_iters,po=po,seed=seed)
    # main(domain_size=2*33,support_size=35,solid_fraction=0.4,h5file="cylinder66.h5",h5group="cylinder66",num_iters=num_iters,po=po,seed=seed)
    # main(domain_size=4*33,support_size=35,solid_fraction=0.4,h5file="cylinder132.h5",h5group="cylinder132",num_iters=num_iters,po=po,seed=seed)
    # main(domain_size=8*33,support_size=35,solid_fraction=0.4,h5file="cylinder264.h5",h5group="cylinder264",bb_correction=bb_correction,num_iters=num_iters,po=4,seed=seed)

    # for n in [32,64,128,256]:
        # build_periodic_domain(num_nodes=n,mesh_type="cartesian",support_size=21)

    # for n in [32,64]:
        # build_periodic_domain(num_nodes=n,mesh_type="scattered",support_size=21,display=True)
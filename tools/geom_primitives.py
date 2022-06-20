import numpy as np
import matplotlib.path as mpath
import matplotlib.pyplot as plt
from numpy.matlib import repmat
import poisson_disk as pd
import h5py

def circle_path(center,radius,step,phi0=0.,normals=False):
    """Create a discretized circle path.

    Args:
        center (tuple): Cartesian coordinates of the circle center.
        radius (float): Radius of the circle
        step (float): Spacing between points on the circle.
        phi0 (float, optional): "Zero" angle of the path.
            If phi0 = 0, the first point is placed at from (xc + R, yc)
        normals (bool): Whether we also want the normal vectors

    Returns:
        If normals is True:
            A tuple of the circle path and normals.
        else:
            A circle path.

    """

    # number of points on half circle
    npoints = 2*int(np.pi*radius/step)

    # create nodes
    center = np.asarray(center)
    nodes = repmat(center,npoints+1,1) # one extra point for closed path

    # angles
    phi1 = np.linspace(phi0,phi0+np.pi,npoints/2,endpoint=False)
    phi2 = np.linspace(phi0+np.pi,phi0+2*np.pi,npoints/2,endpoint=False)

    phi = np.hstack((phi1,phi2))

    # generate points
    nodes[:-1,0] += radius*np.cos(phi)
    nodes[:-1,1] += radius*np.sin(phi)

    # last point of path is equal to first
    nodes[-1,:] = nodes[0,:]

    # path codes
    codes = [mpath.Path.MOVETO] + (npoints-1)*[mpath.Path.LINETO] + [mpath.Path.CLOSEPOLY]

    path = mpath.Path(nodes,codes) # matplotlib path object

    if normals:
        # normal vectors (pointing inwards)
        n = np.empty_like(nodes[:-1,:])
        n[:,0] = -np.cos(phi)
        n[:,1] = -np.sin(phi)
        return path, n
    else:
        return path # matplotlib path object


def box_path(bbox=[0,1,0,1],step=0.1):
    """Create a rectangular box path.

    Args:
        bbox (list): A list containing the coordinates of the bounding box.
            The coordinates should be given as [x1,x2,y1,y2], where
            the point (x1,y1) is the lower left corner, and the point (x2,y2)
            is the upper right corner of the box.
        step (float): The spacing between points on the edges of the box.

    Returns:
        A matplotlib path object representing a discrete box path.
    
    """

    lx = bbox[1] - bbox[0]
    ly = bbox[3] - bbox[2]
    nx = int(lx/step)
    ny = int(ly/step)

    side1 = np.linspace(bbox[0],bbox[1],nx-1,endpoint=False)
    side2 = np.linspace(bbox[2],bbox[3],ny-1,endpoint=False)
    side3 = np.linspace(bbox[1],bbox[0],nx-1,endpoint=False)
    side4 = np.linspace(bbox[3],bbox[2],ny,endpoint=True)

    sides = np.zeros((2*nx + 2*(ny-2)+1,2))
    sides[0:nx-1,0] = side1
    sides[0:nx-1,1] = bbox[2]
    sides[nx-1:nx-1+ny-1,0] = bbox[1] 
    sides[nx-1:nx-1+ny-1,1] = side2
    sides[nx-1+ny-1:nx-1+ny-1+nx-1,0] = side3
    sides[nx-1+ny-1:nx-1+ny-1+nx-1,1] = bbox[3]
    sides[nx-1+ny-1+nx-1:,0] = bbox[0]
    sides[nx-1+ny-1+nx-1:,1] = side4

    codes = len(side1)*[mpath.Path.LINETO]+len(side2)*[mpath.Path.LINETO] \
        +len(side3)*[mpath.Path.LINETO]+len(side4)*[mpath.Path.LINETO]

    codes[0] = mpath.Path.MOVETO
    codes[-1] = mpath.Path.CLOSEPOLY

    return mpath.Path(sides,codes)

def circular_porous_domain(center,radius,dx,seed=None):
    """Creates the geometry for a porous cylinder array.

    Args:
        center (tuple): Cartesian coordinates of the cylinder center.
        radius (float): Radius of the cylinder.
        dx (float): Spacing between points on the edge of the cylinder.
        seed (int): Seed for the random number generator.

    Returns:
        (tuple): A representation of the porous cylinder domain; contains:
            pos: A numpy array of the Cartesian coordinates
            types: A numpy array of integer flags with 1 for interior nodes,
                and -1 for boundary nodes.
            bmap: A numpy array of integer flags with -1 for boundary nodes.
            normals: A numpy array of the normal vectors.

    """

    # Create nodes and normal vectors
    path, normals = circle_path(center=center,radius=radius,step=dx,normals=True)
    bnodes = path.vertices[:-1,:] # matplotlib says vertices should not be accesed this way!

    # Object for Poisson disk sampling
    # Adding a small prefactor to the radius would give tighter points
    sampler = pd.pds(width=1.0,height=1.0,radius=dx,max_shots=40,periodic=True,seed=seed)
    sampler.add_seed_points(bnodes)

    # Generate points
    poisson_points = sampler.rvs()

    # Create mask for points inside of circle path
    mask = path.contains_points(poisson_points,radius=0.1*dx) # not sure what radius is for, but seems to work!
    points = poisson_points[~mask,:]

    pos = np.vstack((bnodes,points))
    btypes = np.empty(bnodes.shape[0],dtype=int)
    btypes[:] = -1
    itypes = np.empty(points.shape[0],dtype=int)
    itypes[:] = 1
    types = np.hstack((btypes,itypes))

    bmap = np.zeros_like(types) # boundary nodes
    bmap[types > 0] = -1        # interior nodes

    print(pos.shape)
    print(types.shape)
    print(bmap.shape)
    print(bmap[bmap==0].shape)
    print(normals.shape)

    # print(types.shape)
    return pos, types, bmap, normals

def periodic_domain(bbox=[0.,1.,0.,1.],dx=0.02,seed=None):
    """Create a rectangular periodic domain.

    Args:
        bbox (list): The bounding box of the periodic domain.
            The coordinates should be given as [x1,x2,y1,y2], where
            the point (x1,y1) is the lower left corner, and the point (x2,y2)
            is the upper right corner of the box.
        dx (float): Minimal spacing between points:
        seed (int): Seed for the random number generator.

    Returns:
        (tuple): A point cloud representation of the periodic domain.
            The tuple contains the Cartesian point coordinates,
            the node types (all interior), and the boundary map (all zeros).

    """

    width = bbox[1] - bbox[0]
    height = bbox[3] - bbox[2]

    sampler = pd.pds(width=width,height=height,radius=dx,max_shots=50,periodic=True,seed=seed)

    points = sampler.rvs()
    types = np.ones(points.shape[0],dtype="i4")
    bmap = np.zeros_like(types)

    points[:,0] += bbox[0]
    points[:,1] += bbox[2]

    print(points.shape)
    print(types.shape)
    print(bmap.shape)
    print(points)
    print(types)
    print(bmap)
    return points,types,bmap


    
def save_to_hdf5(filename,name,pos,types,bmap,normals=None,bbox=None,supports=None):
    """Save a point cloud representation to an HDF5 file.

    Args:
        filename (string): The filename with extension.
        name (string): Group name (kind of like a "root" folder in the file).
        pos (ndarray): Array containing the cartesian point positions.
        types (ndarray): Array of integer flags for node types.
        bmap (ndarray): Array of integer flags for boundary nodes.
        normals (ndarray, optional): Array of the normal vectors for the boundary nodes.
        bbox (list): Bounding box of the point cloud.
        supports (list of lists): The adjacency lists of the points in the point cloud.

    """
    with h5py.File(filename,"w") as f:
        grp = f.create_group(name)
        
        grp.attrs["N"] = pos.shape[0]
        grp.create_dataset("pos",data=pos.T,dtype=np.float64)
        grp.create_dataset("types",data=types,dtype='i4')
        grp.create_dataset("bmap",data=bmap,dtype='i4')

        if normals is not None:
            grp.create_dataset("normals",data=normals.T)
        
        if bbox is not None:
            grp.create_dataset("bbox",data=bbox)

        if supports is not None:
            grp.create_dataset("supports",data=supports.T,dtype='i4')

def save2txt(case,nodes,supports):
    """Save a simple point cloud representation to a set of ascii files.

    The nodes will be saved with the file format ".nodes".
    The supports will be saved with the postfix ".graph".

    Args:
        case (string): A case name. 
        nodes (ndarray): Array containing the cartesian point positions.
        supports (list of lists): The adjacency lists of the points in the point cloud.
    """

    with open(case+".nodes",'w') as file:
        file.write("{}\n".format(len(nodes)))
        for xy in nodes:
            file.write("{} {}\n".format(xy[0],xy[1]))

    nedg = 0
    for row in supports:
        nedg += len(row)

    with open(case+'.graph',mode='w') as file:
        file.write("{} {}\n".format(len(supports),nedg))
        for row in supports:
            file.write(" ".join([str(r) for r in row.tolist()]) + "\n")

def save_support(case,supports):
    """Save a simple point cloud representation to a set of ascii files.

    The nodes will be saved with the file format ".nodes".
    The supports will be saved with the postfix ".graph".

    Args:
        case (string): A case name. 
        supports (list of lists): The adjacency lists of the points in the point cloud.
    """

    nedg = 0
    for row in supports:
        nedg += len(row)

    with open(case+'.graph',mode='w') as file:
        file.write("{} {}\n".format(len(supports),nedg))
        for row in supports:
            file.write(" ".join([str(r) for r in row.tolist()]) + "\n")

def flat2torus(lmbda,phi,bbox):
    """Transformation from periodic plane to 3D torus coordinates.
    
    This function should not be used anymore.
    It is a remnant from a previous attempt to construct periodic 
    supports by performing neighbor searches in three dimensional space.
    It was found that the torus topology skews space, resulting in
    innapropriate stencils when transformed back to two-dimensions.
    The function should not be used anymore.

    Args:
        lmbda (array): List of coordinates in first dimension.
        phi (array): List of coordinates in second dimensions.
        bbox (list): Bounding box of the two-dimensional periodic domain.

    Returns
        x, y, z (tuple): A triplet of coordinates on the surface of a torus.

    """

    # Transform to region [-pi,pi]^2 assuming points are initially in [0,1]^2
    bwidth = bbox[1] - bbox[0]
    bheight = bbox[3] - bbox[2]

    lmbda = lmbda - (0.5*bwidth + bbox[0])
    phi = phi - (0.5*bheight + bbox[2])

    lmbda *= 2*np.pi/bwidth
    phi *= 2*np.pi/bheight

    f = 1./(np.sqrt(2.) - np.cos(phi))

    x3 = f*np.cos(lmbda)
    y3 = f*np.sin(lmbda)
    z3 = f*np.sin(phi)

    return x3, y3, z3

def distances(x,y,supports,bbox):
    """Periodic distance calculations between points.

    Args:
        x (1D array): Array of x-coordinates from the point cloud.
        y (1D array): Array of x-coordinates from the point cloud.
        supports (array): Array of adjacency lists of the point cloud.
        bbox (list): Bounding box of the periodic domain.

    Returns:
        (array): An array distances of the same shape as `supports`.
            We use the bounding box coordinates to account for 
            periodicity in the distance calculation.
    """
    width = bbox[1] - bbox[0]
    height = bbox[3] - bbox[2]

    dist = np.empty(supports.shape)

    for i, supp in enumerate(supports):

        sx = x[supp] - x[supp[0]]
        sy = y[supp] - y[supp[0]]

        lrg = sx > 0.5*width
        smlr = sx <= -0.5*width

        sx[lrg] -= width
        sx[smlr] += width        

        lrg = sy > 0.5*height
        smlr = sy <= -0.5*height

        sy[lrg] -= height
        sy[smlr] += height

        dist[i,:] = np.sqrt(sx**2 + sy**2)

    return dist


def hexagonal_lattice(bbox=[0.,1.,0.,1.],dx=0.02,retstep=False):
    """Returns a point cloud with hexagonal point distribution.

    TODO: An error remains in this function! The spatial step
    size selection mechanism is not guaranteed to work properly.

    Args:
        bbox (list): Bounding box of the periodic domain.
            The aspect ratio of the bounding box should satisfy the
            geometric relationships of equilateral triangles.
        dx (float): The spatial step length between points.
        retstep (bool): Whether to return the spatial step.
            If the dimensions of the bbox and dx do not conform,
            this function attempts to use aslightly different step
            size. 

    Returns:
        A tuple containg the points, types, and boundary map of
        the resulting point cloud. If retstep is True, the actual
        spacing between points is also returned.

    """
    w = bbox[1] - bbox[0] # height
    h = bbox[3] - bbox[2] # width

    from math import ceil

    nx = int(w/dx)
    ny = int(h/(0.5*np.sqrt(3)*dx))
    # if ny % 2 != 0:
        # ny = ny - 1

    # x = np.arange(bbox[0],bbox[1],dx)
    # y = np.arange(bbox[2],bbox[3],0.5*np.sqrt(3)*dx)

    x, dx = np.linspace(bbox[0],bbox[1],nx,endpoint=False,retstep=True)
    print("Hexagonal grid dx = {}".format(dx))
    y = np.linspace(bbox[2],bbox[3],ny,endpoint=False)


    # print(x,y)

    X, Y = np.meshgrid(x,y)
    
    X[1::2,:] += 0.5*dx # shift even rows
    # print(X,Y)

    # print(X.size,Y.size)
    x = X.flatten()
    y = Y.flatten()

    pos = np.stack((x,y),axis=-1)
    # print(pos,pos.shape)
    types = np.ones(pos.shape[0],dtype='i4')
    bmap = -np.ones(pos.shape[0],dtype='i4')

    if retstep:
        return (pos, types, bmap), dx
    else:
        return pos, types, bmap


def cartesian_lattice(bbox=[0.,1.,0.,1.],dx=(0.02,0.02),add_half=False,retstep=False):
    """Returns a Cartesian point cloud distribution.

    Args:
        bbox (list): The bounding box of the periodic domain.
        dx (tuple): The spatial step sizes in x- and y-directions.
        add_half (bool): If the point coordinates are shifted by half of the step size.
            Use this to make the point cloud more "finite-volume" like.
        retstep (bool): Whether to return the actual step size used.
            If the bounding box size is not an integer multiple of the step-size,
            this function will replace the step with a slightly larger one.

    Returns:
        The point cloud representation of the Cartesian lattice.
        If retstep is True, the actual step sizes are also returned.
        
    """

    w = bbox[1] - bbox[0] # height
    h = bbox[3] - bbox[2] # width

    from math import ceil

    ddx, ddy = dx 

    nx = int(w/ddx)
    ny = int(h/ddy)

    x, tx = np.linspace(bbox[0],bbox[1],nx,endpoint=False,retstep=True)
    y, ty = np.linspace(bbox[2],bbox[3],ny,endpoint=False,retstep=True)
    print("Cartesian grid dx, dy = {}, {}".format(tx, ty))

    X, Y = np.meshgrid(x,y)
    
    x = X.flatten()
    y = Y.flatten()

    if add_half:
        x += 0.5*tx
        y += 0.5*ty

    pos = np.stack((x,y),axis=-1)
    types = np.ones(pos.shape[0],dtype='i4')
    bmap = -np.ones(pos.shape[0],dtype='i4')

    if retstep:
        return (pos, types, bmap), (tx,ty)
    else:
        return pos, types, bmap

def test_domain():
    """Test driver for some of the geometric point cloud and I/O primitives.
    """

    print("generating domain")
    p, t, bmap, n = circular_porous_domain(center=(0.5,0.5),radius=0.3,dx=0.01,seed=2019)
    print("domain generation finished")
    bbox = np.asarray([0.,1.,0.,1.])

    save_to_hdf5("geometry.h5","cylinder_array",p,t,bmap,n,bbox)

    b = t == -1
    fig, ax = plt.subplots(nrows=1,ncols=1)
    ax.scatter(p[:,0],p[:,1],s=5,marker='o')
    ax.quiver(p[b,0],p[b,1],n[:,0],n[:,1])
    ax.set_aspect('equal')

    x, y, z = flat2torus(p[:,0],p[:,1],bbox)
    from torus_playground import plot_torus
    fig3, ax3 = plot_torus(x,y,z)

    from scipy.spatial import cKDTree
    
    xyz = np.empty((x.shape[0],3))
    xyz[:,0] = x; xyz[:,1] = y; xyz[:,2] = z

    print("building KDtree")
    kdtree = cKDTree(xyz,leafsize=12)
    print("KDtree is ready")

    x0, y0, z0 = x[234], y[234], z[234]
    d, i = kdtree.query([x0,y0,z0],12)
    print(d,i)

    print("computing nearest neighbours")
    dd, ii = kdtree.query(xyz,12)
    print("nereast neighbours finished")

    save_to_hdf5("geometry_supports.h5","cylinder_array",p,t,bmap,n,bbox,supports=ii)

    ax3.scatter(x[i],y[i],z[i],s=7)
    ax.scatter(p[i,0],p[i,1],marker='x')

    plt.show()

def main():
    """Test driver for some of the geometric point cloud primitives.
    """

    dx = 0.04

    path, n = circle_path(center=(0.5,0.5),radius=0.25,step=dx,normals=True)
    # path = box_path(bbox=[0.25,0.25,0.75,0.75],step=0.05)
    nodes = path.vertices
    print(len(nodes),len(n))

    # Create poisson disk sampling
    max_shots = 40
    Poisson = pd.pds(width=1.0,height=1.0,radius=dx,max_shots=max_shots,periodic=True)

    Poisson.add_seed_points(nodes[:-1,:].tolist())

    points = Poisson.rvs()

    # select points inside path
    mask = path.contains_points(points,radius=dx)

    # join boundary (path) points with interior point
    domain = np.vstack((nodes[:-1],points[~mask,:]))
    # types = np.vstack()

    plt.figure()
    plt.scatter(nodes[:,0],nodes[:,1],marker='o')
    plt.quiver(nodes[:-1,0],nodes[:-1,1],n[:,0],n[:,1])
    plt.scatter(points[mask,0],points[mask,1],c='b')
    plt.scatter(points[~mask,0],points[~mask,1],c='r',marker='x')
    plt.axis('equal')


    # repeat domain for fun
    n = domain.shape[0]
    periodic_domain = repmat(domain,9,1)
    periodic_domain[1:n+1,:] += np.array([1,0])
    periodic_domain[1*n+1:2*n+1,:] += np.array([0,1])
    periodic_domain[2*n+1:3*n+1,:] += np.array([-1,0])
    periodic_domain[3*n+1:4*n+1,:] += np.array([0,-1])
    periodic_domain[4*n+1:5*n+1,:] += np.array([1,1])
    periodic_domain[5*n+1:6*n+1,:] += np.array([-1,1])

    plt.figure()
    plt.scatter(periodic_domain[:,0],periodic_domain[:,1])
    plt.hlines([0,1],-1,2)
    plt.vlines([0,1],-1,2)
    plt.axis('equal')


    plt.show()

if __name__ == '__main__':
    main()
    # test_domain()
#!/usr/bin/env python3

import argparse

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix

from scipy.sparse.csgraph import reverse_cuthill_mckee

try:
    from pymetis import part_graph
except ImportError:
    have_pymetis = False
    print("PyMetis is not available.")
else:
    have_pymetis = True

print("Have_pymetis", have_pymetis)

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



class PeriodicMesh:

    def __init__(self,shape,stepsize,stencil=None):

        self.shape = shape
        self.stepsize = stepsize

    def xy_to_i(self,xy):

        x, y = xy
        return y + self.shape[1]*x

    def i_to_xy(self,i):

        x = i // self.shape[1]
        y = i % self.shape[1]
        return (x, y)

    def draw(self):

        return NotImplementedError

    def _periodic_add(self,p,s):

        x = p[0] + s[0]
        y = p[1] + s[1]

        return (x % self.shape[0], y % self.shape[1])


    def grid_adjacency_csr(self,stencil):
        """Get adjacency structure for a regular periodic grid

        Assume each cell cented has the same stencil, which also means
        each cell/node has the same valency.
        """

        n = self.shape[0]*self.shape[1]

        valency = len(stencil)

        rows = np.empty(n+1,dtype=int)
        cols = np.empty(n*valency,dtype=int)

        rows[0] = 0
        for i in range(0,n):

            rows[i+1] = rows[i] + valency
            xy0 = self.i_to_xy(i)

            for j, s in enumerate(stencil): 

                xy = self._periodic_add(xy0,s)
                cols[rows[i]+j] = self.xy_to_i(xy)

        data = np.ones_like(cols,dtype=int)

        A = csr_matrix((data,cols,rows),shape=(n,n))

        return A

    #def partition(self,method=None):


        # For partitioning regular grids we have four strategies
        #  a) partition along trailing dimension
        #  b) partition along both dimensions
        #  c) partition according to graph of nodal connectivity
        #  d) partition according to graph of elemental connectivity
        #
        # Due to the duality of the grid, methods and c) and d) are 
        # probably equal for regular stencils.



def gridgen():


    parser = argparse.ArgumentParser(description='Cartesian mesh generator in 2D')

    parser.add_argument('-n',nargs='+',type=int, metavar='nd',
        help='Number of cells; if only one argument, the same number of cells is used in each direction')

    parser.add_argument('--stepsize',nargs='+',type=float,
        help='Grid stepsize',default=[1.0])

    parser.add_argument('-f', '--folder', metavar='folder', type=str,
        help='Output folder name')

    parser.add_argument('--cell', action='store_true',
        help='Cell-centered grid')

    parser.add_argument('-p', '--partition', type=int, metavar='nparts', default=0,
        help='Partition the grid into nparts, nparts > 1. If omitted no partitioning is performed.')

    args = parser.parse_args()

    print(args)


    stencil = [(0,0),(1,0),(0,1),(-1,0),(0,-1)] #(1,1),(-1,1),(-1,-1),(1,-1)]

    shape = args.n if len(args.n) > 1 else args.n * 2

    mesh = PeriodicMesh(
        shape=shape,
        stepsize=args.stepsize)

    A = mesh.grid_adjacency_csr(stencil)

    if args.partition:
        ncuts, partition = part_graph(args.partition,xadj=A.indptr,adjncy=A.indices)
        print(partition)

        part_grid = np.reshape(partition,shape)

        xl = np.arange(shape[0]) + 0.5
        yl = np.arange(shape[1]) + 0.5

        x, y = np.meshgrid(xl,yl)

        plt.scatter(x,y,c=partition)
        plt.axis('equal')

        print(A)
        #plt.spy(A)

        # only symmetric, if the stencil is symmetric
        perm = reverse_cuthill_mckee(A,symmetric_mode=True)
        #B = A.tocoo()

        #B.row = perm.take(B.row)
        #B.col = perm.take(B.col)
        #B = B.tocsr()

        #plt.spy(B)

        plt.show()

if __name__ == '__main__':

    gridgen()


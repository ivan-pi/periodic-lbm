import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
from cycler import cycler

import argparse

linestyles = ['solid','dashed','dashdot','dotted']
colors = ['tab:blue','tab:orange','tab:green','tab:red']

class HydroField:

    def __init__(self,filename):

        self.filename = filename

        self._mf = np.load(self.filename)

        if self.is_grid():
            self.rho = self._mf[:,:,0] 
            self.ux  = self._mf[:,:,1] 
            self.uy  = self._mf[:,:,2]
        else:
            self.rho = self._mf[:,0] 
            self.ux  = self._mf[:,1] 
            self.uy  = self._mf[:,2]

    def is_grid(self):
        return self._mf.ndim > 2
    
    @property
    def nfields(self):
        return self._mf.shape[-1]

    @property
    def shape(self):
        return self.rho.shape


def plot_pressure(ax,*files,**kwargs):

    cnt_kwargs = {}

    #cnt_kwargs['colors'] = 'k' if not kwargs['color'] else None
    cnt_kwargs['linewidths'] = kwargs['linewidth'] 

    print(cnt_kwargs)



    for f, ls, clr in zip(files,linestyles,colors):

        hf = HydroField(f)

        if hf.is_grid():

            nx, ny = hf.shape

            x = np.arange(nx) + 0.5
            y = np.arange(ny) + 0.5

            cnt = ax.contour(x,y,hf.rho,linestyles=ls,colors=clr,**cnt_kwargs)

        elif (rho.ndim == 1):

            n, = hf.shape

            raise NotImplementedError


    if kwargs['colorbar']:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        
        plt.colorbar(cnt, cax=cax)


def main():


    parser = argparse.ArgumentParser(description='Vortex post-processor')

    parser.add_argument('files', metavar='file', type=str, nargs='+',
        help='Input file containing vortex fields.')

    parser.add_argument('-o', '--output', metavar='file', type=str,
        help='Output file name')

    parser.add_argument('--figsize', type=float, nargs=2, metavar=('w','h'), default=[4.8,4.8],
        help='Figure size (width, height) in inches, default: (4.8, 4.8)')

    group = parser.add_mutually_exclusive_group()

    group.add_argument('--vbox', type=float, nargs=3, metavar=('x','y','w'),
        help='View box specifying a square region (x-center, y-center, width)')

    group.add_argument('--bbox', type=float, nargs=4, metavar=('l','b', 'w', 'h'),
        help='Bounding box specifying a rectangular region (left, bottom, width, height)')

    parser.add_argument('-n','--levels',type=int, metavar='N',
        help='Number of contour lines')

    parser.add_argument('--lw',type=float, metavar='w',
        help='Contour linewidth, default: 1.5')

    parser.add_argument('--vmin',type=float,
        help='Minimum for colorbar range')

    parser.add_argument('--vmax',type=float,
        help='Maximum for colorbar range')

    parser.add_argument('--caption',type=str,
        help='Figure caption')

    parser.add_argument('--dpi',type=int,default=400,
        help='Dots per inch (dpi), default: 400', )

    parser.add_argument('--color', action='store_true',
        help='Plot colored contours')

    parser.add_argument('--colorbar',action='store_true',
        help='If present, add colorbar to the plot')

    parser.add_argument('--transparent',action='store_true',
        help='Save output with transparent background')

    parser.add_argument('--no-axis',action='store_true')

#
#   parser.add_argument('--nc'

    args = parser.parse_args()

    print(args)

    print("colorbar: {}".format(args.colorbar))

    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=args.figsize)
    
    ax.set_aspect('equal')

    if args.no_axis:
        ax.set_axis_off()

    plot_pressure(ax,*args.files,
        linewidth=args.lw,
        color=args.color,
        colorbar=args.colorbar)

    #
    # Constrain view region
    #
    if args.vbox or args.bbox:

        if args.vbox:
            
            halfwidth = args.vbox[2]/2.

            left  = args.vbox[0] - halfwidth
            right = args.vbox[0] + halfwidth

            bottom = args.vbox[1] - halfwidth
            top    = args.vbox[1] + halfwidth

        elif args.bbox:

            left   = args.bbox[0]
            right  = args.bbox[0] + args.bbox[2]
            
            bottom = args.bbox[1]
            top    = args.bbox[1] + args.bbox[3]

        ax.set_xlim(left, right)
        ax.set_ylim(bottom, top)

    plt.tight_layout()

    #
    # Save figure
    #
    if args.output:

        fig.savefig(
            fname=args.output,
            dpi=args.dpi,
            transparent=args.transparent)

    plt.show()


if __name__ == '__main__':

    main()
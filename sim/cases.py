
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from pathlib import Path
import ctypes

import time
from tqdm import tqdm
import sys

class SimPlugin:


    def __init__(self,path):

        # The underlying C plugin object
        self.ptr = None

        _name = Path(path).stem.removeprefix("lib")
        print("Opening plugin " + _name)

        self.lib = ctypes.CDLL(path)

        # Unpack these for easier referencing
        self.c_init = self.lib["c_"+_name+"_init"]
        self.c_vars = self.lib["c_"+_name+"_vars"]
        self.c_step = self.lib["c_"+_name+"_step"]
        self.c_free = self.lib["c_"+_name+"_free"]

        self.c_norm = self.lib["c_"+_name+"_norm"]
        
        _type = float
        
        self.c_norm.restype = ctypes.c_double
        self.c_norm.argtypes = (ctypes.c_int, ctypes.c_int,
            np.ctypeslib.ndpointer(_type,ndim=2,flags="f_contiguous"),
            np.ctypeslib.ndpointer(_type,ndim=2,flags="f_contiguous"),
            )


        self.c_init.restype = ctypes.c_void_p
        self.c_init.argtypes = (
            ctypes.c_int, 
            ctypes.c_int, 
            ctypes.c_double,
            np.ctypeslib.ndpointer(_type,ndim=2,flags="f_contiguous"),
            np.ctypeslib.ndpointer(_type,ndim=3,flags="f_contiguous"),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_void_p
            ) 

        self.c_step.restype = None
        self.c_step.argtypes = (ctypes.c_void_p, ctypes.c_double)


        self.c_vars.restype = None
        self.c_vars.argtypes = (
            ctypes.c_void_p,
            np.ctypeslib.ndpointer(_type,ndim=2,flags="f_contiguous,writeable"),
            np.ctypeslib.ndpointer(_type,ndim=3,flags="f_contiguous,writeable")
            )

        self.c_free.restype = None
        self.c_free.argtypes = (ctypes.c_void_p,)


    def my_norm(self,u,ua):
        _type = np.float64
        nx, ny = u.shape
        u = np.require(u,_type,['F_CONTIGUOUS'])
        ua = np.require(ua,_type,['F_CONTIGUOUS'])
        return self.c_norm(nx,ny,u,ua)


    def init(self,grid_size,dt,rho,u,sigma=None,params=None):

        self.grid_size = grid_size
        self._dt = dt

        _type = float

        nx, ny = grid_size
        rho = np.require(rho,_type,['F_CONTIGUOUS'])
        u = np.require(u,_type,['F_CONTIGUOUS'])

        _dt = ctypes.c_double(dt)

        if sigma:
            sigma = np.require(sigma, _type, ['F_CONTIGUOUS'])
            sigma_ptr = sigma.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            self.ptr = self.c_init(nx, ny, _dt, rho, u, sigma_p, params)
        else:
            self.ptr = self.c_init(nx, ny, _dt, rho, u, None, params)

    @property
    def dt(self):
        return self._dt

    def vars(self):

        rho = np.empty(self.grid_size,dtype=float,order='F')
        u = np.empty(self.grid_size+(2,),dtype=float,order='F')

        if self.ptr is not None:
            self.c_vars(self.ptr,rho,u)
            return rho, u
        else:
            raise Warning("The plugin has not been initalized. Call .init() first")
            return None, None

    def free(self):

        if self.ptr is not None:
            self.c_free(self.ptr)

    def step(self,omega):

        _omega = ctypes.c_double(omega)
        if self.ptr is not None:
            self.c_step(self.ptr,_omega)

def _vnorm(v):
    return np.linalg.norm(v)
#    return np.sum(v**2)

class TGVStudy:

    # For now we assume we are only studying regular structured grids

    def __init__(self,grid_size,umax,k,nu):
        
        self.umax = umax
        self.nu = nu
        self.k = k

        nx, ny = grid_size
        
        endpoint = False
        x = np.linspace(0,nx,num=nx,endpoint=endpoint)
        y = np.linspace(0,ny,num=ny,endpoint=endpoint)

        x += 0.5
        y += 0.5

        # N.b., indexing arg needed to assure Fortran indexing
        self._xv, self._yv = np.meshgrid(x,y,indexing='ij')

    @property
    def decay_time(self):
        """The vortex decay constant"""
        kx, ky = self.k
        return 1.0/(self.nu*(kx**2 + ky**2))

    @property
    def grid(self):
        return self._xv, self._yv


    def solution(self,t,ret_stress_tensor=False):

        td = self.decay_time
        kx, ky = self.k
        umax = self.umax

        # Aliases
        x = self._xv
        y = self._yv

        u = np.empty(x.shape + (2,),order='F')
        p = np.empty(x.shape,order='F')
      
        u[:,:,0] = -umax*np.sqrt(ky/kx)*np.cos(kx*x)*np.sin(ky*y)*np.exp(-t/td)
        u[:,:,1] =  umax*np.sqrt(kx/ky)*np.sin(kx*x)*np.cos(ky*y)*np.exp(-t/td)
      
        p[:,:]   = -0.25*(umax**2)*((ky/kx)*np.cos(2.0*kx*x) + (kx/ky)*np.cos(2.0*ky*y))*np.exp(-2.0*t/td)

        if ret_stress_tensor:
            S = np.empty(x.shape + (3,),order='F')

            S[:,:,0] = umax*np.sqrt(kx*ky)*np.sin(kx*x)*np.sin(ky*y)*np.exp(-t/td)
            S[:,:,1] = 0.5*umax*(np.sqrt(kx**3/ky) - np.sqrt(ky**3/kx))*np.cos(kx*x)*np.cos(ky*y)*np.exp(-t/td)
            S[:,:,2] = -S[:,:,0]

            return p, u, S

        return p, u

    def norm1(self,t,u):
        """L1-norm of velocity magnitude at time t"""
        _, uref = self.solution(t)
        uerr = u - uref
        magnitude = np.hypot(uerr[:,:,0],uerr[:,:,1])
        return np.max(magnitude)

    def norm2(self,t,u):
        """L2-norm of velocity field at time t"""
        _, uref = self.solution(t)
        mag = np.hypot(u[:,:,0],u[:,:,1])
        magref = np.hypot(uref[:,:,0],uref[:,:,1])
        return _vnorm(magref - mag)/_vnorm(magref)

    def norm2x(self,t,ux):
        """L2-norm of x-velocity at time t"""
        _, uref = self.solution(t)
        return _vnorm(ux - uref[:,:,0])/_vnorm(uref[:,:,0])

    def norm2y(self,t,uy):
        """L2-norm of x-velocity at time t"""
        _, uref = self.solution(t)
        return _vnorm(uy - uref[:,:,1])/_vnorm(uref[:,:,1])


class ShearWaveStudy:

    def __init__(self):
        raise NotImplementedError


def run(sim,omega,nsteps):
    
    start = time.time()
    tm = 0
    
    # Time-stepping loop
    for k in tqdm(range(nsteps)):
        sim.step(omega)
        tm += sim.dt

    elapsed = time.time() - start

    return tm, elapsed


def main(sim):

    from math import ceil


    cs = np.sqrt(1.0/3.0)

    base = 10
    nref = 2

    l1err = np.empty(nref)
    l2err = np.empty(nref)
    l2xerr = np.empty(nref)
    l2yerr = np.empty(nref)

    dt_over_tau = float(sys.argv[1])

    for i in range(nref):

        n = base * 2**i
 
        grid_size = (n, n)
        #umax = 0.01 * cs / (n / base)
        #umax = 0.01  / (n / base)**2
        umax = 0.01  / (n / base)

        k = (2*np.pi/n, 2*np.pi/n)

        Re = 1
        #Re = 100
        nu = umax * n / Re

        tau = nu / cs**2
        
        #cfl = 0.25
        #dt = cfl/np.sqrt(2.0)

        # Corresponds to LBM setting
        #dt = 4 * tau      # Values used by Strzelcyzk
        #dt = 10*tau
        dt = dt_over_tau*tau
        cfl = np.sqrt(2)*dt
        
        omega = dt/(tau + 0.5*dt)

        case = TGVStudy(grid_size,umax=umax,k=k,nu=nu)
        tend = -np.log(0.1) * case.decay_time
        #tend = -np.log(0.5) * case.decay_time
        #tend = -np.log(0.9) * case.decay_time
        nsteps = round(tend / dt)
        #nsteps = 10

        print("CFL    = ", cfl)
        print("Mach   = ", umax/cs)
        print("umax   = ", umax)
        print("nu     = ", nu)
        print("tau    = ", tau)
        print("taulb  = ", 1.0/omega)
        print("omega  = ", omega)
        print("dt     = ", dt)
        print("dt/tau = ", dt/tau)
        print("Re     = ", umax * n / nu)
        print("nsteps = ", nsteps)

        print("Running for {} steps on grid {}".format(nsteps,grid_size))

        p0, u0 = case.solution(t=0.0)
        sim.init(grid_size,dt,p0,u0)

        tend, elapsed = run(sim,omega,nsteps)

        #tend = float(nsteps)
        print("tend = ", tend)

        rho, u = sim.vars()
        _, uref = case.solution(t=tend)

        unorm = np.hypot(u[:,:,0],u[:,:,1])

        print("max u", unorm.max(), (np.hypot(u0[:,:,0],u0[:,:,1])).max())

        l1err[i] = case.norm1(tend,u)
        l2err[i] = case.norm2(tend,u)

        l2xerr[i] = sim.my_norm(u[:,:,0],uref[:,:,0])
        l2xerr[i] = case.norm2x(tend,u[:,:,0])
        l2yerr[i] = case.norm2y(tend,u[:,:,1])

        print(grid_size,l1err[i],l2err[i],elapsed, (n*n*nsteps / elapsed) * 1.0E-6)

        sim.free()

    print("L2: ", l2err)
    print("L2x: ", l2xerr)
    print("L2y: ", l2yerr)

    ux, uy = u[:,:,0], u[:,:,1]
    x, y = case.grid

    urefnorm = np.hypot(uref[:,:,0],uref[:,:,1])
    
    fig, ax = plt.subplots()
    ax.pcolormesh(x,y,unorm)

    ax.quiver(x,y,ux,uy)
    ax.set_aspect('equal')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Taylor-Green Vortex")
    plt.show()

    sz = [base * 2**i for i in range(nref)]

    print("L2, L1")
    for i in range(nref-1):
        print(i, l2err[i]/ l2err[i+1], l1err[i]/l1err[i+1])

    fig, ax = plt.subplots()
    ax.set_xlabel("Nx = Ny")
    ax.set_ylabel("Error")

    fsz = [10,20,40,80,160,320]
    ferr = [
        7.006E-02,
        1.741E-02,
        4.339E-03,
        1.084E-03,
        2.709E-04,
        6.772E-05,
        ]

    #ax.plot(sz,l1err,'s',label=r"$L_1$-error")
    ax.plot(fsz,ferr,label=r"$L_2$-error (nint, nsteps)")
    ax.plot(sz,l2err,'o',label=r"$L_2$-error")
    ax.plot(sz,l2xerr,'v',label=r"$L_2$-x")
    ax.plot(sz,l2yerr,'^',label=r"$L_2$-y")

    szt = np.logspace(1.0,2.5,num=50)
    ax.plot(szt,10.0/szt**2,'--',label=r"$O(1/N^2)$")
    ax.plot(szt,1000./szt**4,'--',label=r"$O(1/N^4)$")
    ax.plot(szt,100000./szt**6,'--',label=r"$O(1/N^6)$")

    ax.set_xscale("log", base=10)
    ax.set_yscale("log", base=10)

    ax.legend()

    plt.show()
    
def _mem9(nx,ny,nhalo=1):
    return (nx+2*nhalo)*(ny+2*nhalo) * 2 * 9 * 8 / (1024**2)

def perfsweep(ax,sim,lo,hi,step=20,tsteps=1000,dt=None,**kwargs):

    grids = list(range(lo,hi,step)) 

    elapsed = np.empty(len(grids),dtype=float)
    mlups = np.empty_like(elapsed)

    cs = np.sqrt(1.0/3.0)

    for i, n in enumerate(grids):

        grid_size = (n, n)
        print("Running grid size {}; est. memory = {} MiB".format(grid_size,_mem9(n,n)))

        umax = 0.01 / (n / lo)
        k = (2*np.pi/n, 2*np.pi/n)
        Re = 100
        nu = umax * n / Re
        tau = nu / cs**2
        if dt is None:
            dt = 1 * tau
        cfl = np.sqrt(2.0)*dt
        omega = dt/(tau + 0.5*dt)
        case = TGVStudy(grid_size,umax=umax,k=k,nu=nu)
        p0, u0 = case.solution(t=0.0)
        sim.init(grid_size,dt,p0,u0)
        del p0, u0

        _, _elapsed = run(sim,omega,tsteps)

        elapsed[i] = _elapsed
        mlups[i] = (n * n / 1000000) * (tsteps/_elapsed)
        sim.free()

    print("grids = ", grids)
    print("mlups = ", mlups)

    ax.plot(grids,mlups,'o-')


def perf():

    matplotlib.use('TkAgg')

    cache = { "L1": 384,"L2": 4*1024, "L3": 16*1024 }

    fig, ax = plt.subplots()

    ax.set_xlabel("N = Nx = Ny")
    ax.set_ylabel("Performance [MLUPS]")
    _bytes_per_cell = 9*2*8

    def n2mem(x):
        return x**2 * _bytes_per_cell / 1024
    def mem2n(x):
        # mem in KiB, 8 bytes per double, square root for equivalent mesh size
        return np.sqrt(x * 1024 / _bytes_per_cell)
    
    nelems = [mem2n(v) for _, v in cache.items()]

    ymin, ymax = ax.get_ylim()
    for nel in nelems:
        ax.axvline(x=nel)

    ax2 = ax.secondary_xaxis("top",transform=(n2mem,mem2n))
    ax2.set_xticks(nelems,labels=cache.keys())

    #sim = SimPlugin("./build/libslbm.so")
    #perfsweep(ax,sim,32,512,step=16,tsteps=500,dt=1.0)
    #perfsweep(ax,sim,512,1024,step=64,tsteps=400,dt=1.0)

    sim = SimPlugin("./build/liblwcuf.so")
    #perfsweep(ax,sim,32,512,step=16,tsteps=500)
    perfsweep(ax,sim,512,1024,step=64,tsteps=400)
    #perfsweep(ax,sim,1024,2048,step=128,tsteps=200)
    #perfsweep(ax,sim,2048,8192,step=256,tsteps=100)

    #sim = SimPlugin("./build/liblw4.so")
    #perfsweep(ax,sim,32,512,step=16,tsteps=500)
    #perfsweep(ax,sim,512,1024,step=64,tsteps=400)

    #sim = SimPlugin("./build/liblw6.so")
    #perfsweep(ax,sim,32,512,step=16,tsteps=400)
    #perfsweep(ax,sim,512,1024,step=64,tsteps=200)
    #perfsweep(ax,sim,1024,2048,step=128,tsteps=100)
    #perfsweep(ax,sim,2048,4096,step=256,tsteps=50)
    #perfsweep(ax,sim,4096,8192,step=256,tsteps=40)
    #perfsweep(ax,sim,8192,9216,step=128,tsteps=40)
    #perfsweep(ax,sim,9216,11264,step=128,tsteps=20)

    ax.set_xscale("log")

    #plt.show()

if __name__ == '__main__':
    

    #sim = SimPlugin("./build/liblwcuf.so")
    #main(sim)

    perf()

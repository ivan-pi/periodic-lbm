
import numpy as np
import matplotlib.pyplot as plt

from importlib import import_module
from ctypes import cdll

"""
The analysis should be encapsulated by a class.

We just drop in our simulation method as a C-compatible function.

The user also gets to initialize any data they may need (say if they want to 
store the fields on the GPU, or add padding, etc.). This data is conserved across calls. 

Essentially we are creating a type of plugin.

We just need to agree on the interface, and what data comes from Python.

Initialization happens here.

The class instance emits the plots we want.

Children of the class implement the various benchmark cases.
Pre-known cases, from the literature are made available for
people to replicate.

They can also expose other kinds of plots (e.g. fitting
the viscosity).

This way, anyone can immediately replicate our study,
just by plugging in their time-stepping method. 

Resulting data can be pickled and shared for comparison 
among researchers.

How can we implement slightly different behavior, with
scattered points, and with padded grids?

"""

def lbm_init_eq()
	"""Standard equilibrium initialization"""
	raise NotImplementedError

def lbm_init_neq(rho,ux,uy,sigma)
	"""Non-equilibrium initialization using stress-tensor"""
	raise NotImplementedError

def _plugin_adaptors(dll):
	"""Load the shared library and returns Python wrappers"""
	raise NotImplementedError

class Simulation:

	def __init__(self,init,step,fin):

		self.init = init
		self.step = step
		self.fin = fin

	@classmethod
	def using_plugin(this,name):
		return this(*_plugin_adaptors(name))

	def initialize(self,case,init_scheme=lbm_init_eq):
		"""Initialize PDF fields for specific case

		In some cases we may like to pass extra arguments, e.g. in case
		of other collision models, hyper-viscosity tunning, performance
		tweaking, etc.
		"""

		# Get macroscopic variables
		self.rho_, self.ux_, self.uy_, _sigma = case(x,y)

		# Get PDFs
		self.pdf = init_scheme(rho,ux,uy,sigma)
		self.ptr = self.init(n,self.pdf)

	def run(nsteps,tau=1.0):
		"""Run a simulation for given number of timesteps"""
		raise NotImplementedError

	@property
	def rho(self):
		return self._rho

	@property
	def ux(self):
		return self._ux

	@property
	def uy(self):
		return self._uy

	@property
	def uy(self):
		return self.pdf


# Q1: Where do we measure the maximum velocity?
# - A) At the point where we expect it to be the maximum
#      according to the known velocity field.
# - B) Using the maximum value of the discrete field
# - C) Using a parabolic fit of the region near the peak
#
# Q2: How do we set the functional form we are fitting?
#     This is specific to each case
#
def viscosity_err(study,relative=True,tmin=0):

	# Exponential decay, governed by viscosity of the liquid
	def f(t,p1,td):
		return p1*exp(-t/td)

	kx, ky = study.wave_vector()
	nu = study.kinematic_viscosity()

	if study.time_data:
	
		time, umax = study.get_time_data()
	
		mask = time >= tmin
		res = optimize.curve_fit(f, time[mask], umax[mask])

		td = res[1]
		nu_meas = 1.0/(td*(kx**2 + ky**2))
		abserr = nu_meas - nu

		if relative:
			return abserr/nu
		else:
			return abserr

	else:
		raise Error("The study does not have time data available.")

class ConvergenceStudy:





if __name__ == '__main__':
	
	# Load the plugin library
	mylbm = cdll.LoadLibrary("mylbm.so")

	# Create the simulation object
	sim = Simulation.using_plugin(mylbm)

	sim.initialize(case,init_scheme="eq")
	sim.run(100)

	study = TGVStudy(sim,
		grid_size=(32,32),
		umax=umax,
		k=(kx,ky),
		nu=nu,
		tmax=tc)

	# For diffusive scaling, the velocity is halved with
	# each refinement stage
	#
	# Initial grid-size = 32
	#
	def spatial_convergence(sim,*args,rsteps=4):

		# TODO: Add viscosity error here

		sz = np.empty(r,dtype=np.int)
		err1 = np.empty(r)
		err2 = np.empty(r)

		umax = args[0]

		for k in range(rsteps):
		
			u = umax / 2**k
			n = 2**(5+k)
			study = TGVStudy(sim,(n,n),u,*args[1:])
			
			sz[k] = n
			err1[k] = study.maxnorm()
			err2[k] = study.norm2()

		return err1, err2


	n, err1, err2 = spatial_convergence(sim,args)
	h = 1.0/(n)

	fig, ax = plt.subplots()

	ax.set_ylabel("error")
	ax.set_xlabel("step size h")

	ax.loglog(h, err1,label="1-norm")
	ax.loglog(h, err2,label="2-norm")


	def temporal_convergence(dt_over_tau,sim,*args):

		err1 = np.empty_like(dt_over_tau)
		err2 = np.empty_like(dt_over_tau)

		for k, v in enumerate(dt_over_tau):

			dt = tau*v

			err1[k] = study.norm1()
			err2[k] = study.norm2()

		return err1, err2

	terr1, terr2 = temporal_convergence()

	fig, ax = plt.subplots()

	ax.set_title("Temporal error")

	ax.set_ylabel("error")
	ax.set_xlabel("dt/tau")

	ax.loglog(dt_over_tau, terr1,label="1-norm")
	ax.loglog(dt_over_tau, terr2,label="2-norm")








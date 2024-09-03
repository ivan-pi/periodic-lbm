module helpers
   use, intrinsic :: iso_c_binding
   implicit none

   interface
      subroutine print_progress(count,max) bind(c)
         import c_size_t
         integer(c_size_t), value :: count, max
      end subroutine
   end interface

end module

program standard_lbm

use kinds, only: wp

use sim_class, only: sim
use sim_fvm_class, only: sim_fvm

use lbm_primitives, only: csqr, lbm_macros, lbm_eqinit, fluid_state
use taylor_green

use helpers, only: print_progress

implicit none

#define DT_OVER_TAU 1

type :: cmd
   integer :: nx, ny, nsteps
#if DT_OVER_TAU
   real(wp) :: dt_over_tau
#endif
end type

type(cmd) :: myargs

type(taylor_green_vortex) :: tgv

class(sim), allocatable :: mysim
real(Wp), allocatable :: rho(:,:), u(:,:,:)

integer :: step, nx, ny, nsteps, write_interval
real(wp) :: omega, dt, nu, tau, nrm, tlbm, umax, kx, ky, tc, tend, ma, cfl

real(kind(1.0d0)) :: start, elapsed

real(wp), parameter :: cs = sqrt(csqr)

myargs = args()

nx = myargs%nx
ny = myargs%ny

! Mach = umax / cs
umax = 0.01_wp * cs / (real(nx,wp) / 64.0_wp)
ma = umax / cs

!umax = ma * cs
! nu = (umax * L) / Re
nu = (umax * real(nx,wp)) / 100._wp
   
!nu = 0.1_wp
!umax = 1.0_wp * nu / real(nx,wp)
!ma = umax / cs

! tau = nu / cs**2
tau = nu / csqr

print *, "Ma = ", ma
print *, "umax = ", umax
print *, "nu   = ", nu
print *, "tau  = ", tau

! Regular LBM
#if DT_OVER_TAU
dt = myargs%dt_over_tau * tau
#else
dt = 5.0*tau
#endif

omega = dt/(tau + 0.5_wp*dt)
cfl = sqrt(2.0_wp)*dt

print *, "tau' = ", 1.0_wp/omega
print *, "omega = ", omega
print *, "dt/tau = ", dt/tau
print *, "CFL = ", cfl

kx = 2*pi/real(nx,wp)
ky = 2*pi/real(ny,wp)

tgv = taylor_green_vortex(k1=kx,k2=ky,umax=umax,nu=nu)
tc = tgv_decay_time(tgv)

call tgv_gnuplot(tgv,"tgv.gp")

! The vortex decays as umax * exp (-t / tc)
! It reaches 10 % of it's initial magnitude at the time

!tend = -tc * log( 0.1_wp )
tend = -tc * log( 0.9_wp )

print *, "tend = ", tend

! Round to nearest integer step
!nsteps = nint(tend / dt)

nsteps = nint(tend / dt)
write_interval = 100

!nsteps = 20000
!write_interval = nsteps + 100
print *, "nsteps = ", nsteps

! Padding calculation
!ldx = nx
!if (mod(nx,16) > 0) ldx = ldx + (16 - mod(nx,16))

allocate(rho(nx,ny))
allocate(u(nx,ny,2))

allocate(sim_fvm :: mysim)

call tgv_ic_field(nx,ny,rho,u)

call mysim%init(nx,ny,dt,rho,u)
call mysim%vars(rho,u)
call output(nx,ny,rho,u(:,:,1),u(:,:,2))

open(22,file="log.txt")

start = wtime()

tlbm = 0
do step = 0, nsteps

      if (mod(step,write_interval) == 0) then
         select type(mysim)
         type is (sim_fvm)
            call mysim%vars(rho,u)
!            call lbm_macros(nx,ny,2,mysim%f1,rho,u(:,:,1),u(:,:,2))
            !call output(nx,ny,lbm%rho,lbm%u,lbm%v, step)
            nrm = norm(u(:,:,1),u(:,:,2),tlbm)
            !print *, tlbm, nrm, maxval(hypot(u(:,:,1),u(:,:,2)))
            write(22,'(3(ES24.17,:,2X))') tlbm, nrm, maxval(hypot(u(:,:,1),u(:,:,2)))
         end select
      end if

      call mysim%step(omega)
      tlbm = tlbm + dt

      call print_progress(int(step,8), int(nsteps,8))
end do

elapsed = wtime() - start
print *, "tend = ", tlbm
print *, "Elapsed = ", elapsed
print *, "MLUPS = ", mlups(nsteps,nx*ny,elapsed)

! Since tlbm was incremented in the last iteration, 
! we need to measure the error using field f2
call mysim%vars(rho,u)
call output(nx,ny,rho,u(:,:,1),u(:,:,2),step)
nrm = norm(u(:,:,1),u(:,:,2),tlbm)

!nrm = norm(lbm%u,lbm%v,tend)

print *, tlbm, nrm, maxval(hypot(u(:,:,1),u(:,:,2)))
write(22,'(3(ES24.17,:,2X))') tlbm, nrm, maxval(hypot(u(:,:,1),u(:,:,2)))
close(22)

print *, nrm

contains

   ! Mega-Lattice Updates per Second
   !
   ! The performance metric is the product of the number of timesteps with the 
   ! number of spatial sites (cells), divided by the elapsed time
   !
   function mlups(nsteps,nsites,elapsed)
      integer, intent(in) :: nsteps,nsites
      real(kind(1.0d0)), intent(in) :: elapsed
      real(kind(1.0d0)) :: mlups
      mlups = (real(nsites,wp) * real(nsteps,wp) * 1.0D-6) / elapsed
   end function

   ! Initial condition callback, as expected by the
   ! equilibrium initialization
   subroutine tgv_ic_field(nx,ny,rho,u)
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: rho(nx,ny), u(nx,ny,2)
      
      integer :: i, j
      real(wp) :: x, y
      type(fluid_state) :: state

      do j = 1, ny
         do i = 1, nx

            x = (i-1) + 0.5_wp
            y = (j-1) + 0.5_wp
            state = tgv_state(tgv,x,y,t=0.0_wp)

            rho = state%p
            u(i,j,1) = state%u
            u(i,j,2) = state%v
         end do
      end do
   end subroutine


   function norm(u,v,t)
      real(wp), intent(in) :: u(:,:), v(:,:), t
      real(wp), allocatable :: ua(:,:), va(:,:)
      real(wp) :: norm

      type(fluid_state) :: state
      real(wp) :: x, y
      integer :: i, j

      allocate(ua,mold=u)
      allocate(va,mold=v)

!      !$omp parallel do collapse(2) if (nx*ny > nomp) &
!      !$omp default(private) shared(nx,ny,ua,va,tgv,t)
      do j = 1, ny
         do i = 1, nx

            x = (i-1) + 0.5_wp
            y = (j-1) + 0.5_wp

            state = tgv_state(tgv,x,y,t)

            ua(i,j) = state%u
            va(i,j) = state%v

         end do
      end do
!      !$omp end parallel do

      norm = norm2(u - ua) / norm2(ua)
      !norm = sqrt(sum((u - ua)**2) / sum(ua**2))
!      norm = norm2(hypot(u-ua,v-va))/norm2(hypot(ua,va))

   end function

   subroutine output(nx,ny,rho,u,v,step)
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: rho(nx,ny), u(nx,ny), v(nx,ny)
      integer, intent(in), optional :: step

      real(wp) :: x, y
      integer :: i, j, wunit
      character(len=60) :: fname

      if (present(step)) then
         write(fname,'("output",I0.7,".dat")') step
      else
         fname = "output.dat"
      end if

      open(newunit=wunit,file=trim(fname),status="unknown",action="write")

      do j = 1, ny
         do i = 1, nx

            x = (i-1) + 0.5_wp
            y = (j-1) + 0.5_wp

            write(wunit,'(5(ES24.17,:,1X))') x, y, rho(i,j), u(i,j), v(i,j)

         end do
      end do

      close(wunit)

   end subroutine


   impure function args()
      type(cmd) :: args

      character(len=24) :: str

      if (command_argument_count() /= 3) then
         write(*,'(A)') "usage: lbm NX NY ratio"
         stop
      end if

      call get_command_argument(1,str)
      read(str,*) args%nx

      call get_command_argument(2,str)
      read(str,*) args%ny

#if DT_OVER_TAU
      call get_command_argument(3,str)
      read(str,*) args%dt_over_tau
#endif

!      call get_command_argument(3,str)
!      read(str,*) args%nsteps

      ! TODO:
      !   - select collision model
      !   - select ddf or normal populations
      !   - initialization (velocity only, pressure, and shear stress)
      !   - output norm
      !   - DUGKS, Simple LBM (tau = 1)
      !   - Thread-safe LBM

   end function


   function wtime()
   !$ use omp_lib, only: omp_get_wtime
      real(kind(1.0d0)) :: wtime
#ifdef _OPENMP 
      wtime = omp_get_wtime()
#else
      call cpu_time(wtime)
#endif
   end function

end program
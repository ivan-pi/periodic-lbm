program standard_lbm

use kinds, only: wp
use lbm_primitives

use taylor_green

implicit none

type :: cmd
   integer :: nx, ny, nsteps
end type

type(cmd) :: myargs

type(taylor_green_vortex) :: tgv

type(lbm_grid(:,:)), allocatable, target :: lbm

integer :: npadded, ldx

integer :: step, nx, ny, nsteps, write_interval
real(wp) :: omega, dt, nu, tau, nrm, tlbm, umax, kx, ky, tc, tend, ma, cfl

real(kind(1.0d0)) :: start, elapsed

real(wp), parameter :: cs = sqrt(csqr)

myargs = args()

nx = myargs%nx
ny = myargs%ny

! Mach = umax / cs
umax = 0.01_wp / (real(nx,wp) / 10.0_wp)
ma = umax / cs

!umax = ma * cs
! nu = (umax * L) / Re
nu = (umax * real(nx,wp)) / 1._wp
   
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
dt = 1.0_wp
omega = 1.0_wp/(tau + 0.5_wp)
cfl = 1.0_wp

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

tend = -tc * log( 0.1_wp )

print *, "tend = ", tend

! Round to nearest integer step
!nsteps = nint(tend / dt)

! Round down
nsteps = nint(tend / dt)
write_interval = 1

!nsteps = 20000
!write_interval = nsteps + 100
print *, "nsteps = ", nsteps

! Padding calculation
!ldx = nx
!if (mod(nx,16) > 0) ldx = ldx + (16 - mod(nx,16))

#if defined(__GFORTRAN__) || defined(__NVCOMPILER)
lbm%nx = nx
lbm%ny = ny
allocate(lbm%f1(0:nx+1,0:ny+1,0:8))
allocate(lbm%f2(0:nx+1,0:ny+1,0:8))
allocate(lbm%rho(nx,ny))
allocate(lbm%u(nx,ny))
allocate(lbm%v(nx,ny))
#else
allocate(lbm_grid(nx,ny) :: lbm)
#endif

call lbm_eqinit(nx,ny,lbm%f1,tgv_ic,nu)

!call lbm_neqinit(nx,ny,lbm%f1,tgv_ic,1.0_wp/omega)

call lbm_macros(nx,ny,lbm%f1,lbm%rho,lbm%u,lbm%v)
call output(nx,ny,lbm%rho,lbm%u,lbm%v)

open(22,file="log.txt")

start = wtime()

tlbm = 0
do step = 1, nsteps
   block
#if defined(__GFORTRAN__) || defined(__NVCOMPILER)
      real(wp), allocatable :: ftmp(:,:,:)

      !call lbm_periodic_bc_pull(nx,ny,lbm%f1)
!      call lbm_collide_and_stream_fused(nx,ny,lbm%f1,lbm%f2,omega)
!      call lbm_periodic_bc_push(nx,ny,lbm%f2)

      call lw_stream(nx,ny,grid%f1,grid%f2,this%dt)
      call lw_collision(nx,ny,grid%f2,omega)
      call lw_bc(nx,ny,grid%f2)

#else
      logical :: even
      even = mod(step,2) == 0
      if (even) then
         call lbm_collide_and_stream_fused(nx,ny,lbm%f1,lbm%f2,omega)
         call lbm_periodic_bc_push(nx,ny,lbm%f2)
      else
         call lbm_collide_and_stream_fused(nx,ny,lbm%f2,lbm%f1,omega)
         call lbm_periodic_bc_push(nx,ny,lbm%f1)
      end if
#endif
      !call lbm_macros(nx,ny,f1,lbm%rho,lbm%u,lbm%v)
      !nrm = norm(lbm%u,lbm%v,tlbm)

      if (mod(step,write_interval) == 0) then
         call lbm_macros(nx,ny,lbm%f1,lbm%rho,lbm%u,lbm%v)
         !call output(nx,ny,lbm%rho,lbm%u,lbm%v, step)
         nrm = norm(lbm%u,lbm%v,tlbm)
         print *, tlbm, nrm, maxval(hypot(lbm%u,lbm%v))
         write(22,'(3(ES24.17,:,2X))') tlbm, nrm, maxval(hypot(lbm%u,lbm%v))
      end if

#if defined(__GFORTRAN__) || defined(__NVCOMPILER)
      call move_alloc(from=lbm%f1,to=ftmp)
      call move_alloc(from=lbm%f2,to=lbm%f1)
      call move_alloc(from=ftmp,to=lbm%f2)
#endif

      tlbm = tlbm + dt

   end block
end do

elapsed = wtime() - start
print *, "tend = ", tlbm
print *, "Elapsed = ", elapsed
print *, "MLUPS = ", mlups(nsteps,nx*ny,elapsed)
print *, "BW (GB/s) ", bandwidth(nsites=nx*ny,nsteps=nsteps,elapsed=elapsed)

! Since tlbm was incremented in the last iteration, 
! we need to measure the error using field f2
call lbm_macros(nx,ny,lbm%f1,lbm%rho,lbm%u,lbm%v)
call output(nx,ny,lbm%rho,lbm%u,lbm%v,step)
nrm = norm(lbm%u,lbm%v,tlbm)

!nrm = norm(lbm%u,lbm%v,tend)

print *, tlbm, nrm, maxval(hypot(lbm%u,lbm%v))
write(22,'(3(ES24.17,:,2X))') tlbm, nrm, maxval(hypot(lbm%u,lbm%v))
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
   pure function tgv_ic(x,y) result(res)
      real(wp), intent(in) :: x, y
      type(fluid_state) :: res

      res = tgv_state(tgv,x,y,t=0.0_wp)

   end function

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

      if (command_argument_count() /= 2) then
         write(*,'(A)') "usage: lbm NX NY"
         stop
      end if

      call get_command_argument(1,str)
      read(str,*) args%nx

      call get_command_argument(2,str)
      read(str,*) args%ny

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

   function bandwidth(nsites,nsteps,elapsed) result(bw)
      integer, parameter :: dp = kind(1.0d0)
      integer, intent(in) :: nsites, nsteps
      real(dp), intent(in) :: elapsed
      real(dp) :: bw 
      
      real(dp), parameter :: bytes_per_gb = 1000.0d0**3

      ! Standard LBM
      !bw = (real(nsites,dp) * 9 * 8 * 2 / bytes_per_gb) / (elapsed / real(nsteps,dp))

      ! LAX WENDROFF
      bw = (real(nsites,dp) * 9 * 8 * 2 * 2 / bytes_per_gb) / (elapsed / real(nsteps,dp))

      ! DUGKS
      !bw = (real(nsites,dp) * 9 * 8 * 2 * 3 / bytes_per_gb) / (elapsed / real(nsteps,dp))

   end function

end program

!> Module for Lax-Wendroff based streaming

! This module has been developed primarily on a RTX 2060 GPU (cc7.5, 12 GB memory)
!
! Relevant resources:
! - https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#compute-capability-7-x
! - https://pde-on-gpu.vaw.ethz.ch/lecture10/
! - https://developer.nvidia.com/blog/using-shared-memory-cuda-cc/
! - https://stackoverflow.com/questions/8011376/when-is-cudas-shared-memory-useful
! - https://developer.nvidia.com/blog/how-access-global-memory-efficiently-cuda-fortran-kernels/
! - https://stackoverflow.com/questions/45532715/optimizing-the-cuda-kernel-with-stencil-pattern
!
! The usage of shared memory for accelerating stencil computations is 
! one of the "standard" usages according to R. Crovella from Nvidia. 

module lbm_lw

use sim_class, only: wp => dp

implicit none
private

public :: lw_collision
public :: lw_stream
public :: lw_bc
public :: lw_vec, lw_vec_fixed

! CUDA kernels
public :: lwcuf_collision
public :: lwcuf_stream

integer, parameter :: cx(0:8) = [0,1,0,-1,0,1,-1,-1,1]
integer, parameter :: cy(0:8) = [0,0,1,0,-1,1,1,-1,-1]

real(wp), parameter :: w0 = 4._wp / 9._wp, &
                       ws = 1._wp / 9._wp, &
                       wd = 1._wp / 36._wp
real(wp), parameter :: w(0:8) = [w0,ws,ws,ws,ws,wd,wd,wd,wd]

real(wp), parameter :: rho0 = 1.0_wp

!integer, parameter :: vxx(0:8) = cx**2
!integer, parameter :: vxy(0:8) = cx*cy
!integer, parameter :: vyy(0:8) = cy**2

type :: lw_vec
   real(wp) :: x(0:8)
   real(wp) :: y(0:8)
   real(wp) :: xx(0:8)
   real(wp) :: xy(0:8)
   real(wp) :: yy(0:8)
end type

contains

   function lw_vec_fixed(dt) result(v)
      real(wp), value :: dt
      type(lw_vec) :: v
      integer :: k
      do k = 1, 8
         v%x(k) = dt*cx(k)
         v%y(k) = dt*cy(k)
         ! dt**2 is carried here
         v%xx(k) = 0.5_wp*v%x(k)*v%x(k)
         v%yy(k) = 0.5_wp*v%y(k)*v%y(k)
         v%xy(k) = v%x(k)*v%y(k)
      end do
   end function

   subroutine lw_macros(nx,ny,fsrc,rho,u,v)
      integer, value :: nx, ny
      real(wp), intent(in) :: fsrc(0:nx+1,0:ny+1,0:8)
      real(wp), intent(out) :: rho(nx,ny), u(nx,ny), v(nx,ny)

      real(wp) :: f(0:8)
      integer :: i, j, k

      ! TODO: ideally, we could somehow overlap the calculation and device
      !       to host copying, for instance by tiling the index space somehow
      !       and launching memcpyAsync

      do concurrent(i=1:nx,j=1:ny)

            ! Gather PDFs
            do k = 0, 8
               f(k) = fsrc(i,j,k)
            end do

            ! density
            rho(i,j) = (((f(5) + f(7)) + (f(6) + f(8))) + &
                      ((f(1) + f(3)) + (f(2) + f(4)))) + f(0)

            rho(i,j) = rho(i,j) + rho0

            ! velocity
            u(i,j) = (((f(5) - f(7)) + (f(8) - f(6))) + (f(1) - f(3))) / rho(i,j)
            v(i,j) = (((f(5) - f(7)) + (f(6) - f(8))) + (f(2) - f(4))) / rho(i,j)

      end do

   end subroutine


   subroutine lw_stream(nx,ny,fsrc,fdst,dt)
      implicit none
      integer, intent(in), value :: nx, ny
      real(wp), intent(in), device :: fsrc(0:nx+1,0:ny+1,0:8)
      real(wp), intent(out), device :: fdst(0:nx+1,0:ny+1,0:8)
      real(wp), intent(in), value :: dt

      ! Streaming variables
      real(wp) :: dfx, dfy, dfxx, dfxy, dfyy

      integer :: i, j, k
      
      real(wp) :: vx(0:8)
      real(wp) :: vy(0:8)
      real(wp) :: vxx(0:8)
      real(wp) :: vxy(0:8)
      real(wp) :: vyy(0:8)

      do k = 1, 8
         vx(k) = dt*cx(k)
         vy(k) = dt*cy(k)
         ! dt**2 is carried here
         vxx(k) = 0.5_wp*vx(k)*vx(k)
         vyy(k) = 0.5_wp*vy(k)*vy(k)
         vxy(k) = vx(k)*vy(k)
      end do

      !$acc data present(fsrc,fdst)

      do concurrent(i=1:nx,j=1:ny)
         fdst(i,j,0) = fsrc(i,j,0)
      end do

      do concurrent(k=1:8,j=1:ny,i=1:nx) &
         local(dfx,dfy,dfxx,dfyy,dfxy)
         
         dfx = 0.5_wp*(fsrc(i+1,j,k) - fsrc(i-1,j,k))
         dfy = 0.5_wp*(fsrc(i,j+1,k) - fsrc(i,j-1,k))

         dfxx = fsrc(i+1,j,k) - 2.0_wp*fsrc(i,j,k) + fsrc(i-1,j,k)
         dfyy = fsrc(i,j+1,k) - 2.0_wp*fsrc(i,j,k) + fsrc(i,j-1,k)
         dfxy = 0.25_wp*(fsrc(i+1,j+1,k) - fsrc(i-1,j+1,k) + fsrc(i-1,j-1,k) - fsrc(i+1,j-1,k))

         fdst(i,j,k) = fsrc(i,j,k) &
            - (vx(k)*dfx + vy(k)*dfy) &
            + (vxx(k)*dfxx + vxy(k)*dfxy + vyy(k)*dfyy)
      end do

      !$acc end data

   end subroutine


   !> This kernel doesn't take care of the rest particles
   !> Those need to be copied separately!
   attributes(global) subroutine lwcuf_stream(nx,ny,fsrc,fdst,vx,vy,vxx,vxy,vyy)
      integer, value :: nx, ny
      real(wp), intent(in) :: fsrc(0:nx+1,0:ny+1,0:8)
      real(wp), intent(out) :: fdst(0:nx+1,0:ny+1,0:8)
      real(wp), intent(in), dimension(0:8) :: vx,vy,vxx,vxy,vyy

      real(wp) :: dfx, dfy, dfxx, dfxy, dfyy
      integer :: i, j, k

#if USE_SHARED
      real(wp), shared :: fs(0:BLOCK_X+1,0:BLOCK_Y+1)
      integer :: is, js
#endif

      i = (blockIdx%x - 1)*blockDim%x + threadIdx%x
      j = (blockIdx%y - 1)*blockDim%y + threadIdx%y
      k = threadIdx%z

#if USE_SHARED
      is = threadIdx%x
      js = threadIdx%y

      ! Make a copy in shared memory

      ! BULK cells
      if (i < (nx+1)+1 .and. j < (ny+1)+1) &
         fs(is,js) = fsrc(i,j,k)

      ! HALO cells
      ! TODO: Fill the layer of HALO cells along the entire edge of the block

      call syncthreads
#endif

      if (i < nx+1 .and. j < ny+1) then
#if USE_SHARED
         dfx = 0.5_wp*(fs(is+1,js) - fs(is-1,js))
         dfy = 0.5_wp*(fs(is,js+1) - fs(is,js-1))
         dfxx = fs(is+1,js) - 2.0_wp*fs(is,js) + fs(is-1,js)
         dfyy = fs(is,js+1) - 2.0_wp*fs(is,js) + fs(is,js-1)
         dfxy = 0.25_wp*(fs(is+1,js+1) - fs(is-1,js+1) + fs(is-1,js-1) - fs(is+1,js-1))
#else
         dfx = 0.5_wp*(fsrc(i+1,j,k) - fsrc(i-1,j,k))
         dfy = 0.5_wp*(fsrc(i,j+1,k) - fsrc(i,j-1,k))
         dfxx = fsrc(i+1,j,k) - 2.0_wp*fsrc(i,j,k) + fsrc(i-1,j,k)
         dfyy = fsrc(i,j+1,k) - 2.0_wp*fsrc(i,j,k) + fsrc(i,j-1,k)
         dfxy = 0.25_wp*(fsrc(i+1,j+1,k) - fsrc(i-1,j+1,k) + fsrc(i-1,j-1,k) - fsrc(i+1,j-1,k))
#endif
         fdst(i,j,k) = fsrc(i,j,k) &
            - (vx(k)*dfx + vy(k)*dfy) &
            + (vxx(k)*dfxx + vxy(k)*dfxy + vyy(k)*dfyy)
      end if

   end subroutine

   attributes(global) subroutine lwcuf_collision(nx,ny,pdf,omega)
      integer, intent(in), value :: nx, ny
      real(wp), intent(inout) :: pdf(0:nx+1,0:ny+1,0:8)
      real(wp), intent(in), value :: omega

      integer :: i, j
      real(wp) :: f(0:8), feq(0:8), rho, ux, uy, uxx, uyy, uxpy, uxmy, indp

      i = (blockIdx%x - 1)*blockDim%x + threadIdx%x
      j = (blockIdx%y - 1)*blockDim%y + threadIdx%y

      if (i > nx .or. j > ny) return

      f(0) = pdf(i,j,0)
      f(1) = pdf(i,j,1)
      f(2) = pdf(i,j,2)
      f(3) = pdf(i,j,3)
      f(4) = pdf(i,j,4)
      f(5) = pdf(i,j,5)
      f(6) = pdf(i,j,6)
      f(7) = pdf(i,j,7)
      f(8) = pdf(i,j,8)

      ! density
      rho = f(0) + (((f(5) + f(7)) + (f(6) + f(8))) + &
             ((f(1) + f(3)) + (f(2) + f(4)))) + rho0

      ! velocity
      ux = (((f(5) - f(7)) + (f(8) - f(6))) + (f(1) - f(3))) / rho
      uy = (((f(5) - f(7)) + (f(6) - f(8))) + (f(2) - f(4))) / rho

      uxx = ux*ux
      uyy = uy*uy
      
      indp = -1.5_wp * (uxx + uyy)
      
      feq(0) = w0*rho*(indp)
      feq(1) = ws*rho*(indp + 3.0_wp*ux + 4.5_wp*uxx)
      feq(2) = ws*rho*(indp + 3.0_wp*uy + 4.5_wp*uyy)
      feq(3) = ws*rho*(indp - 3.0_wp*ux + 4.5_wp*uxx)
      feq(4) = ws*rho*(indp - 3.0_wp*uy + 4.5_wp*uyy)

      uxpy = ux + uy
      feq(5) = wd*rho*(indp + 3.0_wp*uxpy + 4.5_wp*uxpy*uxpy)
      feq(7) = wd*rho*(indp - 3.0_wp*uxpy + 4.5_wp*uxpy*uxpy)
      
      uxmy = ux - uy
      feq(6) = wd*rho*(indp - 3.0_wp*uxmy + 4.5_wp*uxmy*uxmy)
      feq(8) = wd*rho*(indp + 3.0_wp*uxmy + 4.5_wp*uxmy*uxmy)

      feq = feq + w*(rho - rho0)

      pdf(i,j,0) = f(0) + omega*(feq(0) - f(0))
      pdf(i,j,1) = f(1) + omega*(feq(1) - f(1))
      pdf(i,j,2) = f(2) + omega*(feq(2) - f(2))
      pdf(i,j,3) = f(3) + omega*(feq(3) - f(3))
      pdf(i,j,4) = f(4) + omega*(feq(4) - f(4))
      pdf(i,j,5) = f(5) + omega*(feq(5) - f(5))
      pdf(i,j,6) = f(6) + omega*(feq(6) - f(6))
      pdf(i,j,7) = f(7) + omega*(feq(7) - f(7))
      pdf(i,j,8) = f(8) + omega*(feq(8) - f(8))

   end subroutine

   subroutine lw_collision(nx,ny,pdf,omega)
      implicit none
      integer, intent(in), value :: nx, ny
      real(wp), intent(inout), device :: pdf(0:nx+1,0:ny+1,0:8)
      real(wp), intent(in), value :: omega

      ! Collision variables
      real(wp) :: rho, ux, uy, feq(0:8), f(0:8)
      real(wp) :: indp, uxx, uyy, uxpy, uxmy

      integer :: i, j

      !$acc data present(pdf)

      do concurrent(j=1:ny,i=1:nx) &
         local(rho,ux,uy,uxx,uyy,indp,uxpy,uxmy,feq,f)

            ! Gather PDFs
            f(0) = pdf(i,j,0)
            f(1) = pdf(i,j,1)
            f(2) = pdf(i,j,2)
            f(3) = pdf(i,j,3)
            f(4) = pdf(i,j,4)
            f(5) = pdf(i,j,5)
            f(6) = pdf(i,j,6)
            f(7) = pdf(i,j,7)
            f(8) = pdf(i,j,8)

            ! density
            rho = f(0) + (((f(5) + f(7)) + (f(6) + f(8))) + &
                   ((f(1) + f(3)) + (f(2) + f(4)))) + rho0

            ! velocity
            ux = (((f(5) - f(7)) + (f(8) - f(6))) + (f(1) - f(3))) / rho
            uy = (((f(5) - f(7)) + (f(6) - f(8))) + (f(2) - f(4))) / rho

            uxx = ux*ux
            uyy = uy*uy
            
            indp = -1.5_wp * (uxx + uyy)
            
            feq(0) = w0*rho*(indp)
            feq(1) = ws*rho*(indp + 3.0_wp*ux + 4.5_wp*uxx)
            feq(2) = ws*rho*(indp + 3.0_wp*uy + 4.5_wp*uyy)
            feq(3) = ws*rho*(indp - 3.0_wp*ux + 4.5_wp*uxx)
            feq(4) = ws*rho*(indp - 3.0_wp*uy + 4.5_wp*uyy)

            uxpy = ux + uy
            feq(5) = wd*rho*(indp + 3.0_wp*uxpy + 4.5_wp*uxpy*uxpy)
            feq(7) = wd*rho*(indp - 3.0_wp*uxpy + 4.5_wp*uxpy*uxpy)
            
            uxmy = ux - uy
            feq(6) = wd*rho*(indp - 3.0_wp*uxmy + 4.5_wp*uxmy*uxmy)
            feq(8) = wd*rho*(indp + 3.0_wp*uxmy + 4.5_wp*uxmy*uxmy)

            feq = feq + w*(rho - rho0)

            pdf(i,j,0) = f(0) + omega*(feq(0) - f(0))
            pdf(i,j,1) = f(1) + omega*(feq(1) - f(1))
            pdf(i,j,2) = f(2) + omega*(feq(2) - f(2))
            pdf(i,j,3) = f(3) + omega*(feq(3) - f(3))
            pdf(i,j,4) = f(4) + omega*(feq(4) - f(4))
            pdf(i,j,5) = f(5) + omega*(feq(5) - f(5))
            pdf(i,j,6) = f(6) + omega*(feq(6) - f(6))
            pdf(i,j,7) = f(7) + omega*(feq(7) - f(7))
            pdf(i,j,8) = f(8) + omega*(feq(8) - f(8))

      end do

      !$acc end data

   end subroutine


   subroutine lw_bc(nx,ny,fdst)
      implicit none
      integer, intent(in), value :: nx, ny
      real(wp), intent(inout), device :: fdst(0:nx+1,0:ny+1,0:8)

      integer :: i, j

      ! In the Lax-Wendroff method, we need to couple all of the 
      ! PDFs and not just the incoming ones!

      !$acc data present(fdst)

      ! EAST and WEST HALO
      do concurrent(j=1:ny)
         fdst(0,j,:) = fdst(nx,j,:)
         fdst(nx+1,j,:) = fdst(1,j,:)
      end do

      ! NORTH and SOUTH HALO
      do concurrent(i=0:nx+1)
         fdst(i,0,:) = fdst(i,ny,:)
         fdst(i,ny+1,:) = fdst(i,1,:)
      end do

      !$acc end data

   end subroutine

end module



!> Module implementing standard LBM method
module sim_lw_class

   use cudafor   
   use sim_class, only: sim, dp
   
   use lbm_primitives, only: lbm_eqinit, lbm_macros
   use lbm_lw, only: lw_collision, lw_stream, lw_bc, lw_vec, lw_vec_fixed

   implicit none

   private
   public :: sim_lw

   type, extends(sim) :: sim_lw
      ! nx and ny are part of the parent
      real(dp), allocatable, device :: f1(:,:,:), f2(:,:,:)
      real(dp) :: dt
      type(lw_vec) :: v
 
      ! Internal variables for performance measurement
      integer(8) :: clock_rate
      real(8) :: stm = 0, ctm = 0
      integer :: nsteps = 0
   contains
      procedure :: init => lw_init
      procedure :: step => lw_step
      procedure :: vars => lw_vars
   end type

contains

   subroutine lw_init(this,nx,ny,dt,rho,u,sigma)
      class(sim_lw), intent(out) :: this
      integer, intent(in) :: nx, ny
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: rho(nx,ny), u(nx,ny,2)    ! rho, (u, v)
      real(dp), intent(in), optional :: sigma(nx,ny,3)  ! (sxx, sxy, syy)

      real(dp), allocatable :: fh(:,:,:)

      ! Store the time-step
      this%dt = dt

      ! Store grid-sizes
      this%nx = nx
      this%ny = ny

      this%v = lw_vec_fixed(dt)

      ! TODO: padding for optimal memory access
      ! - Consider using also routines like CudaMalloc3D

      allocate(this%f1(0:nx+1,0:ny+1,0:8))
      allocate(this%f2(0:nx+1,0:ny+1,0:8))

      ! Temporary host buffer
      allocate(fh(0:nx+1,0:ny+1,0:8))

      ! Initialize the PDFs using the provided fields
      if (present(sigma)) then
         error stop "NotImplementedError"
      else
         ! equilibrium initialization
         call lbm_eqinit(nx,ny,1,fh,rho,u(:,:,1),u(:,:,2))        
      end if

      ! Copy from host to device
      this%f1 = fh

      print *, "maxval(fh) = ", maxval(fh)

      ! Fill the HALO layers
      call lw_bc(nx,ny,this%f1)

      CALL SYSTEM_CLOCK(count_rate=this%clock_rate)

   end subroutine

   subroutine lw_step(this,omega)

      !@cuf use cudafor, only: dim3, cudaStreamSynchronize
      class(sim_lw), intent(inout) :: this
      real(dp), intent(in) :: omega

      real(dp), allocatable, device :: ftmp(:,:,:)

      real(dp), device, dimension(0:8) :: vx, vy, vxx, vxy, vyy
      integer :: istat, i, j
      integer(8) :: tstart, tend
      type(dim3) :: grid, tBlock

#define USE_CUDA_KERNELS 0
      
! ------------------------
      call SYSTEM_CLOCK(tstart)

#if USE_CUDA_KERNELS
      vx = this%v%x
      vy = this%v%y
      vxx = this%v%xx
      vxy = this%v%xy
      vyy = this%v%yy

      tBlock = dim3(64,4,4)
      grid = dim3(ceiling(real(this%nx)/tBlock%x),ceiling(real(this%ny)/tBlock%y), 2)

      do concurrent(i=1:this%nx,j=1:this%ny)
         this%f2(i,j,0) = this%f1(i,j,0)
      end do

      call lwcuf_stream<<<grid,tBlock>>>(this%nx,this%ny,this%f1,this%f2,&
        vx,vy,vxx,vxy,vyy)
      istat = cudaStreamSynchronize()
#else
      call lw_stream(this%nx,this%ny,this%f1,this%f2,this%dt)
#endif 
      call SYSTEM_CLOCK(tend)
      this%stm = this%stm + real(tend - tstart,8)/real(this%clock_rate,8)

      call SYSTEM_CLOCK(tstart)

#if USE_CUDA_KERNELS
      tBlock = dim3(64,4,1)
      grid = dim3(ceiling(real(this%nx)/tBlock%x),ceiling(real(this%ny)/tBlock%y), 1)
      call lwcuf_collision<<<grid,tBlock>>>(this%nx,this%ny,this%f2,omega)
      istat = cudaStreamSynchronize()
#else
      call lw_collision(this%nx,this%ny,this%f2,omega)
#endif
      call SYSTEM_CLOCK(tend)
      this%ctm = this%ctm + real(tend - tstart,8)/real(this%clock_rate,8)

      call lw_bc(this%nx,this%ny,this%f2)

      call move_alloc(from=this%f1,to=ftmp)
      call move_alloc(from=this%f2,to=this%f1)
      call move_alloc(from=ftmp,to=this%f2)

      this%nsteps = this%nsteps + 1

   end subroutine

   subroutine lw_vars(this,rho,u)
      class(sim_lw), intent(in) :: this
      real(dp), intent(out), contiguous :: rho(:,:), u(:,:,:)

      real(dp), allocatable :: fh(:,:,:)
      
      ! Temporary buffer to preserve intent
      allocate(fh(0:this%nx+1,0:this%ny+1,0:8))

      ! COPY DEVICE to HOST
      fh = this%f1

      call lbm_macros(this%nx,this%ny,1,fh,rho,u(:,:,1),u(:,:,2))
   
   end subroutine

end module


!
! PLUGIN INTERFACE
!
!void *siminit(
!   int nx, int ny, double dt,
!   double *rho, double *u, double *sigma, void *params);
function c_lwcuf_init(nx,ny,dt,rho,u,sigma,params) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_int, c_loc
   use sim_lw_class
   implicit none
   integer(c_int), value :: nx, ny
   real(c_double), value :: dt
   real(c_double), intent(in) :: rho(nx,ny), u(nx,ny,2)
   real(c_double), intent(in), optional :: sigma(nx,ny,3)
   type(c_ptr), value :: params
   type(c_ptr) :: c_lwcuf_init

   type(sim_lw), pointer :: sim_ptr 
   
   print *, "Hello in INIT wrapper"

   allocate(sim_ptr)
   call sim_ptr%init(nx,ny,dt,rho,u,sigma)
   c_lwcuf_init = c_loc(sim_ptr)

   print *, "Goodbye in INIT wrapper"

end function

!void simstep(void *sim, double omega);
subroutine c_lwcuf_step(ptr,omega) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_f_pointer
   use sim_lw_class
   type(c_ptr), value :: ptr
   real(c_double), value :: omega

   type(sim_lw), pointer :: sim_ptr
   
   call c_f_pointer(ptr, sim_ptr)
   call sim_ptr%step(omega)

end subroutine

!void simvars(void *sim, double *rho, double *u);
subroutine c_lwcuf_vars(ptr,rho,u) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_f_pointer
   use sim_lw_class
   type(c_ptr), value :: ptr
   real(c_double), intent(out), target :: rho(*), u(*)

   real(c_double), pointer, contiguous :: rho_(:,:), u_(:,:,:)
   type(sim_lw), pointer :: sim_ptr 

   call c_f_pointer(ptr, sim_ptr)
   
   associate(nx => sim_ptr%nx, ny => sim_ptr%ny)

   rho_(1:nx,1:ny) => rho(1:nx*ny)
   u_(1:nx,1:ny,1:2) => u(1:nx*ny*2)
   call sim_ptr%vars(rho_,u_)

   end associate

end subroutine


!void simfree(void *sim);
subroutine c_lwcuf_free(ptr) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_associated
   use sim_lw_class
   implicit none
   type(c_ptr), value :: ptr
   
   !TODO: verify if this can be polymorphic
   type(sim_lw), pointer :: simp 
   
   if(c_associated(ptr)) then
      call c_f_pointer(ptr, simp)
      
      print *, "Streaming BW (GB/s): ", ((sizeof(simp%f1) + sizeof(simp%f2))/1.0D9)/(simp%stm/simp%nsteps) 
      print *, "Collision BW (GB/s): ", ((sizeof(simp%f1) + sizeof(simp%f2))/1.0D9)/(simp%ctm/simp%nsteps) 

      DEALLOCATE(simp)
   end if

end subroutine


function c_lwcuf_norm(nx,ny,u,ua) bind(c)
   use, intrinsic :: iso_c_binding, only: c_int, c_double
   use sim_lw_class
   integer(c_int), value :: nx, ny
   real(c_double), intent(in) :: u(nx,ny), ua(nx,ny)
   real(c_double) :: c_lwcuf_norm
   c_lwcuf_norm = norm2(u - ua) / norm2(ua)
end function


module lbm_fvm

use sim_class, only: wp => dp

implicit none
private

!public :: fvm_step
public :: fvm_collision
public :: fvm_predict_hc
public :: fvm_correct_hc

public :: fvm_bc

real(wp), parameter :: cx(0:8) = [0,1,0,-1,0,1,-1,-1,1]
real(wp), parameter :: cy(0:8) = [0,0,1,0,-1,1,1,-1,-1]

real(wp), parameter :: w0 = 4._wp / 9._wp, &
                       ws = 1._wp / 9._wp, &
                       wd = 1._wp / 36._wp
real(wp), parameter :: w(0:8) = [w0,ws,ws,ws,ws,wd,wd,wd,wd]

real(wp), parameter :: rho0 = 1.0_wp

      ! regular derivative parameters
      real(wp), parameter :: p18 = 0.125_wp, &
                             p34 = 0.75_wp, &
                             p38 = 0.375_wp

contains

   pure function equilibrium(rho,ux,uy) result(feq)
      !$acc routine seq
      real(wp), intent(in) :: rho, ux, uy
      real(wp) :: feq(0:8)

      real(wp) :: uxx, uyy,  uxpy, uxmy
      real(wp) :: indp

      uxx = ux*ux
      uyy = uy*uy

      !
      ! DDF Shifting (for better accuracy)
      !
      
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

   end function


   ! HEUN-central PREDICTOR
   subroutine fvm_predict_hc(nx,ny,fsrc,fdst,dt)
      implicit none
      integer, intent(in), value :: nx, ny
      real(wp), intent(in)  :: fsrc(-1:nx+2,-1:ny+2,0:8)    ! \hat{g}
      real(wp), intent(inout) :: fdst(-1:nx+2,-1:ny+2,0:8)   ! g^*
      real(wp), intent(in), value :: dt

      ! Streaming variables
      integer :: i, j, k

      real(wp) :: flux
      

      ! Zero-population
      !$acc kernels
      fdst(:,:,0) = fsrc(:,:,0)
      !$acc end kernels

      !$acc parallel loop collapse(3) private(flux)
      !$omp parallel do collapse(3) private(flux)
      do k = 1, 8
      do j = 1, ny
         do i = 1, nx

            flux = cx(k)*(0.5_wp*(fsrc(i+1,j,k) + fsrc(i,j,k)) - 0.5_wp*(fsrc(i,j,k) + fsrc(i-1,j,k))) &
                +  cy(k)*(0.5_wp*(fsrc(i,j+1,k) + fsrc(i,j,k)) - 0.5_wp*(fsrc(i,j,k) + fsrc(i,j-1,k)))

            fdst(i,j,k) = fdst(i,j,k) - dt*flux

         end do
      end do
      end do

   end subroutine fvm_predict_hc

   ! HEUN-central CORRECTOR
   subroutine fvm_correct_hc(nx,ny,fp,fsrc,fdst,dt)
      implicit none
      integer, intent(in), value :: nx, ny
      real(wp), intent(in) :: fp(-1:nx+2,-1:ny+2,0:8)
      real(wp), intent(in) :: fsrc(-1:nx+2,-1:ny+2,0:8)
      real(wp), intent(inout) :: fdst(-1:nx+2,-1:ny+2,0:8)
      real(wp), intent(in), value :: dt

      ! Streaming variables
      integer :: i, j, k
      real(wp) :: flux, fluxp

      !$acc parallel loop collapse(3) private(flux,fluxp)
      !$omp parallel do collapse(3) private(flux,fluxp)
      do k = 1, 8
      do j = 1, ny
         do i = 1, nx

            flux = cx(k)*(0.5_wp*(fsrc(i+1,j,k) + fsrc(i,j,k)) - 0.5_wp*(fsrc(i,j,k) + fsrc(i-1,j,k))) &
                +  cy(k)*(0.5_wp*(fsrc(i,j+1,k) + fsrc(i,j,k)) - 0.5_wp*(fsrc(i,j,k) + fsrc(i,j-1,k)))

            fluxp = cx(k)*(0.5_wp*(fp(i+1,j,k) + fp(i,j,k)) - 0.5_wp*(fp(i,j,k) + fp(i-1,j,k))) &
                 +  cy(k)*(0.5_wp*(fp(i,j+1,k) + fp(i,j,k)) - 0.5_wp*(fp(i,j,k) + fp(i,j-1,k)))


            fdst(i,j,k) = fdst(i,j,k) - dt*0.5_wp*(flux + fluxp)

         end do
      end do
      end do

   end subroutine fvm_correct_hc

   subroutine fvm_collision(nx,ny,pdf,omega)
      implicit none
      integer, intent(in), value :: nx, ny
      real(wp), intent(inout) :: pdf(-1:nx+2,-1:ny+2,0:8)
      real(wp), intent(in), value :: omega

      ! Collision variables
      real(wp) :: rho, ux, uy, fneq(0:8), f(0:8)
      
      integer :: i, j, k

      !
      ! N.b. the HALO layer is collided too!
      !

      !$acc parallel loop gang collapse(2) private(rho,ux,uy,fneq,f)
      !$omp parallel do collapse(2) private(f,fneq,rho,ux,uy)
      do j = -1, ny+2
         do i = -1, nx+2

            ! Gather PDFs
            !$acc loop seq
            do k = 0, 8
               f(k) = pdf(i,j,k)
            end do

            ! density
            rho = f(0) + (((f(5) + f(7)) + (f(6) + f(8))) + &
                   ((f(1) + f(3)) + (f(2) + f(4)))) + rho0

            ! velocity
            ux = (((f(5) - f(7)) + (f(8) - f(6))) + (f(1) - f(3))) / rho
            uy = (((f(5) - f(7)) + (f(6) - f(8))) + (f(2) - f(4))) / rho

            fneq = equilibrium(rho, ux, uy) - f

            ! Scatter PDFs
            pdf(i,j,0) = f(0) + omega*fneq(0)
            pdf(i,j,1) = f(1) + omega*fneq(1)
            pdf(i,j,2) = f(2) + omega*fneq(2)
            pdf(i,j,3) = f(3) + omega*fneq(3)
            pdf(i,j,4) = f(4) + omega*fneq(4)
            pdf(i,j,5) = f(5) + omega*fneq(5)
            pdf(i,j,6) = f(6) + omega*fneq(6)
            pdf(i,j,7) = f(7) + omega*fneq(7)
            pdf(i,j,8) = f(8) + omega*fneq(8)

         end do
      end do

   end subroutine


   subroutine fvm_bc(nx,ny,fdst)
      implicit none
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: fdst(-1:nx+2,-1:ny+2,0:8)

      ! In the Lax-Wendroff method, we need to couple all of the 
      ! PDFs and not just the incoming ones!

      !$acc kernels

      ! SOUTH HALO
      fdst(1:nx,   0,:) = fdst(1:nx,ny  ,:)
      fdst(1:nx,  -1,:) = fdst(1:nx,ny-1,:)

      ! NORTH HALO
      fdst(1:nx,ny+1,:) = fdst(1:nx,1,:)
      fdst(1:nx,ny+2,:) = fdst(1:nx,2,:)

      ! WEST HALO
      fdst( 0,-1:ny+2,:) = fdst(nx  ,-1:ny+2,:)
      fdst(-1,-1:ny+2,:) = fdst(nx-1,-1:ny+2,:)

      ! EAST HALO
      fdst(nx+1,-1:ny+2,:) = fdst( 1,-1:ny+2,:)
      fdst(nx+2,-1:ny+2,:) = fdst( 2,-1:ny+2,:)

      !$acc end kernels

   end subroutine

end module


#define HEUN_QUICK 0

!> Module implementing standard LBM method
module sim_fvm_class
   
   use sim_class, only: sim, dp
   
   use lbm_primitives, only: lbm_eqinit, lbm_macros


   use lbm_fvm, only: fvm_collision, fvm_bc, &
#if HEUN_QUICK
      fvm_predict => fvm_predict_quick, &
      fvm_correct => fvm_correct_quick
#else
      fvm_predict => fvm_predict_hc, &
      fvm_correct => fvm_correct_hc
#endif

   implicit none

   type, extends(sim) :: sim_fvm
      real(dp), allocatable :: f1(:,:,:), f2(:,:,:)
      real(dp), allocatable :: fc(:,:,:)
      real(dp) :: dt
   contains
      procedure :: init => fvm_init
      procedure :: step => fvm_step
      procedure :: vars => fvm_vars
   end type

contains

   subroutine fvm_init(this,nx,ny,dt,rho,u,sigma)
      class(sim_fvm), intent(out) :: this
      integer, intent(in) :: nx, ny
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: rho(nx,ny), u(nx,ny,2)    ! rho, (u, v)
      real(dp), intent(in), optional :: sigma(nx,ny,3)  ! (sxx, sxy, syy)

      ! Store the time-step
      this%dt = dt

      ! Store grid size
      this%nx = nx
      this%ny = ny

      ! Allocate memory buffers
      allocate(this%f1(-1:nx+2,-1:ny+2,0:8))
      allocate(this%f2(-1:nx+2,-1:ny+2,0:8))
      allocate(this%fc(-1:nx+2,-1:ny+2,0:8))


      ! Equilibrium initialization
      call lbm_eqinit(nx,ny,2,this%f1,rho,u(:,:,1),u(:,:,2))

      ! TODO: BC!!!
      call fvm_bc(nx,ny,this%f1)

   end subroutine

   subroutine fvm_step(this,omega)
      class(sim_fvm), intent(inout) :: this
      real(dp), intent(in) :: omega

      real(dp) :: tau_g, omega_g

      ! omega := dt/(tau + dt/2)
      ! tau_g = tau/dt + 1/2

      tau_g = 1.0_dp/omega

      ! 1. Collision             g -> \hat{g}
      ! 2. Prediction           \hat{g} -> g^*, \hat{R}^t
      ! 3. Correction           \hat{g} + 

      ! Copy original populations
      this%fc = this%f1
      call fvm_collision(this%nx,this%ny,this%fc,omega) ! g -> \hat{g}

      ! --------------------------
      this%f2 = this%fc
      call fvm_collision(this%nx,this%ny,this%f1,0.5_dp*omega)

      call fvm_predict(this%nx,this%ny,fsrc=this%f1,fdst=this%f2,dt=this%dt)  ! g -= dt*(R^{t}/2), g -> g^*
      call fvm_bc(this%nx,this%ny,this%f2)

      ! f2 now contains f*

      ! -------------------
      call fvm_collision(this%nx,this%ny,this%f2,0.5_dp*omega)
      call fvm_correct(this%nx,this%ny,fp=this%f2,fsrc=this%f1,fdst=this%fc,dt=this%dt)  ! g -= dt*(R^{t+h})
      call fvm_bc(this%nx,this%ny,this%fc)

      swap: block
         real(dp), allocatable :: ftmp(:,:,:)
         call move_alloc(from=this%fc,to=ftmp)
         call move_alloc(from=this%f1,to=this%fc)
         call move_alloc(from=ftmp,to=this%f1)
      end block swap

   end subroutine

   subroutine fvm_vars(this,rho,u)
      class(sim_fvm), intent(in) :: this
      real(dp), intent(out), contiguous :: rho(:,:), u(:,:,:)

      !$acc update self(this%f1)
      call lbm_macros(this%nx,this%ny,2,this%f1,rho,u(:,:,1),u(:,:,2))

   end subroutine

end module




!
! PLUGIN INTERFACE
!
!void *siminit(
!   int nx, int ny, double dt,
!   double *rho, double *u, double *sigma, void *params);
function c_fvm_init(nx,ny,dt,rho,u,sigma,params) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_int, c_loc
   use sim_fvm_class
   implicit none
   integer(c_int), value :: nx, ny
   real(c_double), value :: dt
   real(c_double), intent(in) :: rho(nx,ny), u(nx,ny,2)
   real(c_double), intent(in), optional :: sigma(nx,ny,3)
   type(c_ptr), value :: params
   type(c_ptr) :: c_fvm_init

   type(sim_fvm), pointer :: sim_ptr 
   
   print *, "Hello in INIT wrapper"

   allocate(sim_ptr)
   call sim_ptr%init(nx,ny,dt,rho,u,sigma)
   c_fvm_init = c_loc(sim_ptr)

   print *, "Goodbye in INIT wrapper"

end function

!void simstep(void *sim, double omega);
subroutine c_fvm_step(ptr,omega) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_f_pointer
   use sim_fvm_class
   type(c_ptr), value :: ptr
   real(c_double), value :: omega

   type(sim_fvm), pointer :: sim_ptr
   
   call c_f_pointer(ptr, sim_ptr)
   call sim_ptr%step(omega)

end subroutine

!void simvars(void *sim, double *rho, double *u);
subroutine c_fvm_vars(ptr,rho,u) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_f_pointer
   use sim_fvm_class
   type(c_ptr), value :: ptr
   real(c_double), intent(out), target :: rho(*), u(*)

   real(c_double), pointer, contiguous :: rho_(:,:), u_(:,:,:)
   type(sim_fvm), pointer :: sim_ptr 

   call c_f_pointer(ptr, sim_ptr)
   
   associate(nx => sim_ptr%nx, ny => sim_ptr%ny)

   rho_(1:nx,1:ny) => rho(1:nx*ny)
   u_(1:nx,1:ny,1:2) => u(1:nx*ny*2)
   call sim_ptr%vars(rho_,u_)

   end associate

end subroutine


!void simfree(void *sim);
subroutine c_fvm_free(ptr) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_associated
   use sim_fvm_class, only: sim_fvm
   implicit none
   type(c_ptr), value :: ptr
   
   !TODO: verify if this can be polymorphic
   type(sim_fvm), pointer :: simp 
   
   if(c_associated(ptr)) then
      call c_f_pointer(ptr, simp)
      DEALLOCATE(simp)
   end if

end subroutine


function c_fvm_norm(nx,ny,u,ua) bind(c)
   use, intrinsic :: iso_c_binding, only: c_int, c_double
   use sim_fvm_class
   integer(c_int), value :: nx, ny
   real(c_double), intent(in) :: u(nx,ny), ua(nx,ny)
   real(c_double) :: c_fvm_norm
   c_fvm_norm = norm2(u - ua) / norm2(ua)
end function

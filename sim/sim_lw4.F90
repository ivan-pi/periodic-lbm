module lbm_lw4

use sim_class, only: wp => dp

implicit none
private

public :: lw4_collision
public :: lw4_stream
public :: lw4_bc, nhalo

integer, parameter :: cx(0:8) = [0,1,0,-1,0,1,-1,-1,1]
integer, parameter :: cy(0:8) = [0,0,1,0,-1,1,1,-1,-1]

real(wp), parameter :: w0 = 4._wp / 9._wp, &
                       ws = 1._wp / 9._wp, &
                       wd = 1._wp / 36._wp
real(wp), parameter :: w(0:8) = [w0,ws,ws,ws,ws,wd,wd,wd,wd]

real(wp), parameter :: rho0 = 1.0_wp

integer, parameter :: nhalo = 2

contains

   subroutine lw4_stream(nx,ny,fsrc,fdst,dt)
      implicit none
      integer, intent(in), value :: nx, ny
      real(wp), intent(in)  :: fsrc(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8), dt
      real(wp), intent(out) :: fdst(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8)

      ! Streaming variables
      real(wp) :: dfx, dfy, dfxx, dfxy, dfyy
      real(wp) :: fc, fn, fs, fw, fe, fne, fnw, fse, fsw
      real(wp) :: fnn, fss, fww, fee, fne2, fnw2, fse2, fsw2

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


      !
      ! Streaming
      !
      fdst(1:nx,1:ny,0) = fsrc(1:nx,1:ny,0)


      ! TODO: investigate use of nontemporal clause on fdst

      !$omp parallel do simd collapse(3) schedule(simd: static) &
      !$omp default(private) shared(nx,ny,fsrc,fdst,vx,vy,vxx,vxy,vyy)
      do k = 1, 8
         do j = 1, ny
            do i = 1, nx
               !
               ! Gather PDFs from neighbouring nodes
               !

               fc  = fsrc(i  , j, k)
               
               fe  = fsrc(i+1, j, k)
               fw  = fsrc(i-1, j, k)

               fee  = fsrc(i+2, j, k)
               fww  = fsrc(i-2, j, k)

               fn  = fsrc(i, j+1, k)
               fs  = fsrc(i, j-1, k)
               
               fnn  = fsrc(i, j+2, k)
               fss  = fsrc(i, j-2, k)

               fne = fsrc(i+1, j+1, k)
               fnw = fsrc(i-1, j+1, k)
               fsw = fsrc(i-1, j-1, k)
               fse = fsrc(i+1, j-1, k)

               fne2 = fsrc(i+2, j+2, k)
               fnw2 = fsrc(i-2, j+2, k)
               fsw2 = fsrc(i-2, j-2, k)
               fse2 = fsrc(i+2, j-2, k)

               dfx = (1.0_wp/12.0_wp)*(fww - fee) + (2.0_wp/3.0_wp)*(fe - fw)
               dfy = (1.0_wp/12.0_wp)*(fss - fnn) + (2.0_wp/3.0_wp)*(fn - fs)

               dfxx = -(1.0_wp/12.0_wp)*fww + (4.0_wp/3.0_wp)*fw - (5.0_wp/2.0_wp)*fc  + (4.0_wp/3.0_wp)*fe - (1.0_wp/12.0_wp)*fee
               dfyy = -(1.0_wp/12.0_wp)*fss + (4.0_wp/3.0_wp)*fs - (5.0_wp/2.0_wp)*fc  + (4.0_wp/3.0_wp)*fn - (1.0_wp/12.0_wp)*fnn

               dfxy = (1.0_wp/3.0_wp)*(fne - fnw + fsw - fse) - (1.0_wp/48.0_wp)*(fne2 - fnw2 + fsw2 - fse2)

               !dfxx = fe - 2.0_wp*fc + fw
               !dfyy = fn - 2.0_wp*fc + fs
               !dfxy = 0.25_wp*(fne - fse - fnw + fsw)

               fdst(i,j,k) = fc - vx(k)*dfx - vy(k)*dfy + (vxx(k)*dfxx + vxy(k)*dfxy + vyy(k)*dfyy)

            end do
         end do
      end do

   end subroutine

   subroutine lw4_collision(nx,ny,pdf,omega)
      implicit none
      integer, intent(in), value :: nx, ny
      real(wp), intent(inout) :: pdf(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8)
      real(wp), intent(in) :: omega

      ! Collision variables
      real(wp) :: rho, irho, ux, uy, feq(0:8), f(0:8)
      real(wp) :: indp, uxx, uyy, uxpy, uxmy
      real(wp) :: omegabar

      integer :: i, j

      omegabar = 1.0_wp - omega
      !
      ! Collision
      !

      !$omp parallel do simd collapse(2) schedule(simd: static) &
      !$omp default(private) shared(nx,ny,pdf) firstprivate(omega,omegabar)
      do j = 1, ny
         do i = 1, nx

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

            irho = 1.0_wp/rho

            ! velocity
            ux = (((f(5) - f(7)) + (f(8) - f(6))) + (f(1) - f(3))) * irho
            uy = (((f(5) - f(7)) + (f(6) - f(8))) + (f(2) - f(4))) * irho

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
      end do

   end subroutine


   subroutine lw4_bc(nx,ny,fdst)
      implicit none
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: fdst(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8)

      ! In the Lax-Wendroff method, we need to couple all of the 
      ! PDFs and not just the incoming ones!

#if 1
      ! Copy top and bottom layers
      fdst(1:nx,-1:0,:) = fdst(1:nx,ny-1:ny,:)
      fdst(1:nx,ny+1:ny+2,:) = fdst(1:nx,1:2,:)

      ! Copy left and right layers
      fdst(-1:0,:,:) = fdst(nx-1:nx,-1:ny+2,:)
      fdst(nx+1:nx+2,:,:) = fdst(1:2,:,:)
#else
      ! SOUTH HALO
      fdst(1:nx,0,:) = fdst(1:nx,ny,:)

      ! NORTH HALO
      fdst(1:nx,ny+1,:) = fdst(1:nx,1,:)

      ! WEST HALO
      fdst(0,1:ny,:) = fdst(nx,1:ny,:)

      ! EAST HALO
      fdst(nx+1,1:ny,:) = fdst(1,1:ny,:)

      ! Second layer
      fdst(0:nx+1,  -1,:) = fdst(0:nx+1,ny-1,:)
      fdst(0:nx+1,ny+2,:) = fdst(0:nx+1, 2,:)
      fdst(  -1,0:ny+1,:) = fdst(nx-1,0:ny+1,:)
      fdst(nx+2,0:ny+1,:) = fdst( 2,0:ny+1,:)

      !
      ! CORNERS
      !
       fdst( 0, 0,:) = fdst(nx  , ny  , :) 
       fdst(-1,-1,:) = fdst(nx-1, ny-1, :)

       fdst(nx+1,ny+1,:) = fdst(1, 1,:) 
       fdst(nx+2,ny+2,:) = fdst(2, 2,:) 

       fdst( 0, ny+1, :) = fdst(nx  ,1, :) 
       fdst(-1, ny+2, :) = fdst(nx-1,2, :)

       fdst(nx+1, 0, :) = fdst(1,ny  , :)
       fdst(nx+2,-1, :) = fdst(2,ny-1, :)
#endif
   end subroutine

end module



!> Module implementing standard LBM method
module sim_lw4_class
   
   use sim_class, only: sim, dp
   
   use lbm_primitives, only: &
      lbm_grid, lbm_eqinit, lbm_macros

   use lbm_lw4, only: lw4_collision, lw4_stream, lw4_bc, nhalo
   implicit none

   private
   public :: sim_lw4

   type, extends(sim) :: sim_lw4
#if defined(__GFORTRAN__) || defined(__NVCOMPILER)
      type(lbm_grid) :: grid
#else
      type(lbm_grid(:,:)), allocatable :: grid
#endif
      integer :: flip = 0
      real(dp) :: dt
   contains
      procedure :: init => lw4_init
      procedure :: step => lw4_step
      procedure :: vars => lw4_vars
   end type

contains

   subroutine lw4_init(this,nx,ny,dt,rho,u,sigma)
      class(sim_lw4), intent(out) :: this
      integer, intent(in) :: nx, ny
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: rho(nx,ny), u(nx,ny,2)    ! rho, (u, v)
      real(dp), intent(in), optional :: sigma(nx,ny,3)  ! (sxx, sxy, syy)

      ! Store the time-step
      this%dt = dt

      ! PDTs for the win!
#if defined(__GFORTRAN__) || defined(__NVCOMPILER)
   print *, "Hello in INIT"
      this%grid%nx = nx
      this%grid%ny = ny
      allocate(this%grid%f1(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8))
      allocate(this%grid%f2(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8))
#else
      ! PDTs for the win!
      allocate(lbm_grid(nx,ny,halo=nhalo) :: this%grid)
#endif
      this%flip = 0
!      call move_alloc(from=grid,to=this%grid)


      ! Initialize the PDFs using the provided fields
      if (present(sigma)) then
         ! non-equilibrium initialization
!         call lbm_neqinit(nx,ny,this%grid%f1, &
!            rho,u(:,:,1),u(:,:,2), &
!            sigma(:,:,1),sigma(:,:,2),sigma(:,:,3))
         error stop "NotImplementedError"
      else
         ! equilibrium initialization
         call lbm_eqinit(nx,ny,nhalo,this%grid%f1,rho,u(:,:,1),u(:,:,2))
      end if

      ! Fill the HALO layers
      call lw4_bc(nx,ny,this%grid%f1)

      !this%grid%f2 = this%grid%f1

   end subroutine

   subroutine lw4_step(this,omega)
      class(sim_lw4), intent(inout) :: this
      real(dp), intent(in) :: omega

      associate( &
         nx => this%grid%nx, ny => this%grid%ny, &
         grid => this%grid)

#if defined(__GFORTRAN__) || defined(__NVCOMPILER)
            call lw4_stream(nx,ny,grid%f1,grid%f2,this%dt)
            call lw4_collision(nx,ny,grid%f2,omega)
            call lw4_bc(nx,ny,grid%f2)
            swap: block
               real(dp), allocatable :: ftmp(:,:,:)
               call move_alloc(from=grid%f1,to=ftmp)
               call move_alloc(from=grid%f2,to=grid%f1)
               call move_alloc(from=ftmp,to=grid%f2)
            end block swap
#else
         select case(this%flip)
         case(0)
            call lw4_stream(nx,ny,grid%f1,grid%f2,this%dt)
            call lw4_collision(nx,ny,grid%f2,omega)
            call lw4_bc(nx,ny,grid%f2)
         case(1)
            call lw4_stream(nx,ny,grid%f2,grid%f1,this%dt)
            call lw4_collision(nx,ny,grid%f1,omega)
            call lw4_bc(nx,ny,grid%f1)
         end select
         this%flip = 1 - this%flip
#endif

      end associate

   end subroutine

   subroutine lw4_vars(this,rho,u)
      class(sim_lw4), intent(in) :: this
      real(dp), intent(out), contiguous :: rho(:,:), u(:,:,:)
      associate(grid => this%grid)
#if defined(__GFORTRAN__) || defined(__NVCOMPILER)
         call lbm_macros(grid%nx,grid%ny,nhalo,grid%f1,rho,u(:,:,1),u(:,:,2))
#else
         select case(this%flip)
         case(0)
            call lbm_macros(grid%nx,grid%ny,nhalo,grid%f1,rho,u(:,:,1),u(:,:,2))
         case(1)
            call lbm_macros(grid%nx,grid%ny,nhalo,grid%f2,rho,u(:,:,1),u(:,:,2))
         end select
#endif
      end associate
   end subroutine


end module


!
! PLUGIN INTERFACE
!
!void *siminit(
!   int nx, int ny, double dt,
!   double *rho, double *u, double *sigma, void *params);
function c_lw4_init(nx,ny,dt,rho,u,sigma,params) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_int, c_loc
   use sim_lw4_class
   implicit none (type, external)
   integer(c_int), value :: nx, ny
   real(c_double), value :: dt
   real(c_double), intent(in) :: rho(nx,ny), u(nx,ny,2)
   real(c_double), intent(in), optional :: sigma(nx,ny,3)
   type(c_ptr), value :: params
   type(c_ptr) :: c_lw4_init

   type(sim_lw4), pointer :: sim_ptr 
   
   print *, "Hello in INIT wrapper"

   allocate(sim_ptr)
   call sim_ptr%init(nx,ny,dt,rho,u,sigma)
   c_lw4_init = c_loc(sim_ptr)

   print *, "Goodbye in INIT wrapper"

end function

!void simstep(void *sim, double omega);
subroutine c_lw4_step(ptr,omega) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_f_pointer
   use sim_lw4_class
   type(c_ptr), value :: ptr
   real(c_double), value :: omega

   type(sim_lw4), pointer :: sim_ptr
   
   call c_f_pointer(ptr, sim_ptr)
   call sim_ptr%step(omega)

end subroutine

!void simvars(void *sim, double *rho, double *u);
subroutine c_lw4_vars(ptr,rho,u) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_f_pointer
   use sim_lw4_class
   type(c_ptr), value :: ptr
   real(c_double), intent(out), target :: rho(*), u(*)

   real(c_double), pointer, contiguous :: rho_(:,:), u_(:,:,:)
   type(sim_lw4), pointer :: sim_ptr 

   call c_f_pointer(ptr, sim_ptr)
   
   associate(nx => sim_ptr%grid%nx, ny => sim_ptr%grid%ny)

   rho_(1:nx,1:ny) => rho(1:nx*ny)
   u_(1:nx,1:ny,1:2) => u(1:nx*ny*2)
   call sim_ptr%vars(rho_,u_)

   end associate

end subroutine


!void simfree(void *sim);
subroutine c_lw4_free(ptr) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_associated
   use sim_lw4_class
   implicit none
   type(c_ptr), value :: ptr
   
   !TODO: verify if this can be polymorphic
   type(sim_lw4), pointer :: simp 
   
   if(c_associated(ptr)) then
      call c_f_pointer(ptr, simp)
      DEALLOCATE(simp)
   end if

end subroutine


function c_lw4_norm(nx,ny,u,ua) bind(c)
   use, intrinsic :: iso_c_binding, only: c_int, c_double
   use sim_lw4_class
   integer(c_int), value :: nx, ny
   real(c_double), intent(in) :: u(nx,ny), ua(nx,ny)
   real(c_double) :: c_lw4_norm
   c_lw4_norm = norm2(u - ua) / norm2(ua)
end function


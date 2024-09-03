module lbm_lw4cuf

use sim_class, only: wp => dp

implicit none
private

public :: lw4cuf_collision
public :: lw4cuf_stream
public :: lw4cuf_bc, nhalo

integer, parameter :: cx(0:8) = [0,1,0,-1,0,1,-1,-1,1]
integer, parameter :: cy(0:8) = [0,0,1,0,-1,1,1,-1,-1]

real(wp), parameter :: w0 = 4._wp / 9._wp, &
                       ws = 1._wp / 9._wp, &
                       wd = 1._wp / 36._wp
real(wp), parameter :: w(0:8) = [w0,ws,ws,ws,ws,wd,wd,wd,wd]

real(wp), parameter :: rho0 = 1.0_wp

integer, parameter :: nhalo = 2

contains

   subroutine lw4cuf_stream(nx,ny,fsrc,fdst,dt)
      implicit none
      integer, intent(in), value :: nx, ny
      real(wp), intent(in), device  :: fsrc(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8)
      real(wp), intent(out), device :: fdst(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8)
      real(wp), intent(in), value :: dt

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

      do concurrent(i=1:nx,j=1:ny,k=1:8) &
         local(fc,fe,fw,fee,fww,fn,fs,fnn,fss,&
            dfx,dfy,dfxx,dfxy,dfyy)

         fc  = fsrc(i  , j, k)
         
         fe  = fsrc(i+1, j, k)
         fw  = fsrc(i-1, j, k)

         fee  = fsrc(i+2, j, k)
         fww  = fsrc(i-2, j, k)

         fn  = fsrc(i, j+1, k)
         fs  = fsrc(i, j-1, k)
         
         fnn  = fsrc(i, j+2, k)
         fss  = fsrc(i, j-2, k)

         dfx = (2.0_wp/3.0_wp)*(fsrc(i+1,j,k) - fsrc(i-1,j,k)) - (1.0_wp/12.0_wp)*(fsrc(i+2,j,k) - fsrc(i-2,j,k))
         dfy = (2.0_wp/3.0_wp)*(fsrc(i,j+1,k) - fsrc(i,j-1,k)) - (1.0_wp/12.0_wp)*(fsrc(i,j+2,k) - fsrc(i,j-2,k))

         dfxx = -(1.0_wp/12.0_wp)*fww + (4.0_wp/3.0_wp)*fw - (5.0_wp/2.0_wp)*fc  + (4.0_wp/3.0_wp)*fe - (1.0_wp/12.0_wp)*fee
         dfyy = -(1.0_wp/12.0_wp)*fss + (4.0_wp/3.0_wp)*fs - (5.0_wp/2.0_wp)*fc  + (4.0_wp/3.0_wp)*fn - (1.0_wp/12.0_wp)*fnn

         ! dfxy = (1.0_wp/3.0_wp)*(fsrc(i+1,j+1,k) - fsrc(i-1,j+1,k) + fsrc(i-1,j-1,k) - fsrc(i+1,j-1,k)) &
         !      - (1.0_wp/48.0_wp)*(fsrc(i+2,j+2,k) - fsrc(i-2,j+2,k) + fsrc(i-2,j-2,k) - fsrc(i+2,j-2,k))
         
         dfxy = (5.0_wp/12.0_wp)*(fsrc(i+1,j+1,k) - fsrc(i-1,j+1,k) + fsrc(i-1,j-1,k) - fsrc(i+1,j-1,k)) &
               - (1.0_wp/24.0_wp)*(fsrc(i+1,j+2,k) - fsrc(i-2,j+1,k) + fsrc(i-1,j-2,k) - fsrc(i+2,j-1,k)) &
               - (1.0_wp/24.0_wp)*(fsrc(i+2,j+1,k) - fsrc(i-1,j+2,k) + fsrc(i-2,j-1,k) - fsrc(i+1,j-2,k))

         ! dfxy = (4.0_wp/9.0_wp)*(fsrc(i+1,j+1,k) - fsrc(i-1,j+1,k) + fsrc(i-1,j-1,k) - fsrc(i+1,j-1,k)) &
         !      - (1.0_wp/18.0_wp)*(fsrc(i+1,j+2,k) - fsrc(i-2,j+1,k) + fsrc(i-1,j-2,k) - fsrc(i+2,j-1,k)) &
         !      - (1.0_wp/18.0_wp)*(fsrc(i+2,j+1,k) - fsrc(i-1,j+2,k) + fsrc(i-2,j-1,k) - fsrc(i+1,j-2,k)) &
         !      + (1.0_wp/144.0_wp)*(fsrc(i+2,j+2,k) - fsrc(i-2,j+2,k) + fsrc(i-2,j-2,k) - fsrc(i+2,j-2,k))

         fdst(i,j,k) = fsrc(i,j,k) &
            - (vx(k)*dfx + vy(k)*dfy) &
            + (vxx(k)*dfxx + vxy(k)*dfxy + vyy(k)*dfyy)
      end do

   end subroutine

   subroutine lw4cuf_collision(nx,ny,pdf,omega)
      implicit none
      integer, intent(in), value :: nx, ny
      real(wp), intent(inout), device :: pdf(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8)
      real(wp), intent(in), value :: omega

      ! Collision variables
      real(wp) :: rho, irho, ux, uy, feq(0:8), f(0:8)
      real(wp) :: indp, uxx, uyy, uxpy, uxmy

      integer :: i, j


      do concurrent(i=1:nx,j=1:ny) local(rho,ux,uy,feq,f,indp,uxx,uyy,uxpy,uxmy)

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
            ux = (((f(5) - f(7)) + (f(8) - f(6))) + (f(1) - f(3))) /rho
            uy = (((f(5) - f(7)) + (f(6) - f(8))) + (f(2) - f(4))) /rho

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

   end subroutine


   subroutine lw4cuf_bc(nx,ny,fdst)
      implicit none
      integer, intent(in), value :: nx, ny
      real(wp), intent(inout), device :: fdst(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8)

      integer :: i, j
      ! In the Lax-Wendroff method, we need to couple all of the 
      ! PDFs and not just the incoming ones!

      ! Copy top and bottom layers
      do concurrent(i=1:nx)
         fdst(i,-1,:) = fdst(i,ny-1,:)
         fdst(i, 0,:) = fdst(i,ny   ,:)
         fdst(i,ny+1,:) = fdst(i,1,:)
         fdst(i,ny+2,:) = fdst(i,2,:)
      end do

      ! Copy left and right layers
      do concurrent(j=-1:ny+2)
         fdst(-1,j,:) = fdst(nx-1,j,:)
         fdst( 0,j,:) = fdst(nx  ,j,:)
         fdst(nx+1,j,:) = fdst(1,j,:)
         fdst(nx+2,j,:) = fdst(2,j,:)
      end do

   end subroutine

end module


!> Lax-Wendroff LBM with 4th-order FDM streaming
module sim_lw4cuf_class
   
   use sim_class, only: sim, dp
   
   use lbm_primitives, only: lbm_eqinit, lbm_macros
   use lbm_lw4cuf, only: lw4cuf_collision, lw4cuf_stream, lw4cuf_bc, nhalo
   implicit none

   private
   public :: sim_lw4cuf

   type, extends(sim) :: sim_lw4cuf
      ! nx and ny are part of the parent
      real(dp), allocatable, device :: f1(:,:,:), f2(:,:,:)
      real(dp) :: dt
      integer(8) :: clock_rate
      real(8) :: stm = 0, ctm = 0
      integer :: nsteps = 0
   contains
      procedure :: init => lw4cuf_init
      procedure :: step => lw4cuf_step
      procedure :: vars => lw4cuf_vars
   end type

contains

   subroutine lw4cuf_init(this,nx,ny,dt,rho,u,sigma)
      class(sim_lw4cuf), intent(out) :: this
      integer, intent(in) :: nx, ny
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: rho(nx,ny), u(nx,ny,2)    ! rho, (u, v)
      real(dp), intent(in), optional :: sigma(nx,ny,3)  ! (sxx, sxy, syy)

      real(dp), allocatable :: fh(:,:,:)

      print *, "initialization started"

      ! Store the time-step
      this%dt = dt

      this%nx = nx
      this%ny = ny
      allocate(this%f1(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8))
      allocate(this%f2(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8))

      ! Temporary host buffer
      allocate(fh(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8))

      ! Equilibrium Initialization
      call lbm_eqinit(nx,ny,2,fh,rho,u(:,:,1),u(:,:,2))

      ! Copy from host to device
      this%f1 = fh

      ! Fill the HALO layers
      call lw4cuf_bc(nx,ny,this%f1)

      CALL SYSTEM_CLOCK(count_rate=this%clock_rate)

      print *, "initialization done"

   end subroutine

   subroutine lw4cuf_step(this,omega)
      class(sim_lw4cuf), intent(inout) :: this
      real(dp), intent(in) :: omega
      
      integer(8) :: tstart, tend
      real(dp), device, allocatable :: ftmp(:,:,:)

      call SYSTEM_CLOCK(tstart)

      call lw4cuf_stream(this%nx,this%ny,this%f1,this%f2,this%dt)
      
      call SYSTEM_CLOCK(tend)
      this%stm = this%stm + real(tend - tstart,8)/real(this%clock_rate,8)


      call SYSTEM_CLOCK(tstart)

      call lw4cuf_collision(this%nx,this%ny,this%f2,omega)

      call SYSTEM_CLOCK(tend)
      this%ctm = this%ctm + real(tend - tstart,8)/real(this%clock_rate,8)

      call lw4cuf_bc(this%nx,this%ny,this%f2)


      ! SWAP
      call move_alloc(from=this%f1,to=ftmp)
      call move_alloc(from=this%f2,to=this%f1)
      call move_alloc(from=ftmp,to=this%f2)

      this%nsteps = this%nsteps + 1

   end subroutine

   subroutine lw4cuf_vars(this,rho,u)
      class(sim_lw4cuf), intent(in) :: this
      real(dp), intent(out), contiguous :: rho(:,:), u(:,:,:)

      real(dp), allocatable :: fh(:,:,:)
      
      ! Temporary buffer to preserve intent
      allocate(fh(0:this%nx+1,0:this%ny+1,0:8))

      ! COPY DEVICE to HOST
      fh = this%f1

      call lbm_macros(this%nx,this%ny,nhalo,fh,rho,u(:,:,1),u(:,:,2))

   end subroutine


end module


!
! PLUGIN INTERFACE
!
!void *siminit(
!   int nx, int ny, double dt,
!   double *rho, double *u, double *sigma, void *params);
function c_lw4cuf_init(nx,ny,dt,rho,u,sigma,params) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_int, c_loc
   use sim_lw4cuf_class
   implicit none
   integer(c_int), value :: nx, ny
   real(c_double), value :: dt
   real(c_double), intent(in) :: rho(nx,ny), u(nx,ny,2)
   real(c_double), intent(in), optional :: sigma(nx,ny,3)
   type(c_ptr), value :: params
   type(c_ptr) :: c_lw4cuf_init

   type(sim_lw4cuf), pointer :: sim_ptr 
   
   print *, "Hello in INIT wrapper"

   allocate(sim_ptr)
   call sim_ptr%init(nx,ny,dt,rho,u,sigma)
   c_lw4cuf_init = c_loc(sim_ptr)

   print *, "Goodbye in INIT wrapper"

end function

!void simstep(void *sim, double omega);
subroutine c_lw4cuf_step(ptr,omega) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_f_pointer
   use sim_lw4cuf_class
   type(c_ptr), value :: ptr
   real(c_double), value :: omega

   type(sim_lw4cuf), pointer :: sim_ptr
   
   call c_f_pointer(ptr, sim_ptr)
   call sim_ptr%step(omega)

end subroutine

!void simvars(void *sim, double *rho, double *u);
subroutine c_lw4cuf_vars(ptr,rho,u) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_f_pointer
   use sim_lw4cuf_class
   type(c_ptr), value :: ptr
   real(c_double), intent(out), target :: rho(*), u(*)

   real(c_double), pointer, contiguous :: rho_(:,:), u_(:,:,:)
   type(sim_lw4cuf), pointer :: sim_ptr 

   call c_f_pointer(ptr, sim_ptr)
   
   associate(nx => sim_ptr%nx, ny => sim_ptr%ny)

   rho_(1:nx,1:ny) => rho(1:nx*ny)
   u_(1:nx,1:ny,1:2) => u(1:nx*ny*2)
   call sim_ptr%vars(rho_,u_)

   end associate

end subroutine


!void simfree(void *sim);
subroutine c_lw4cuf_free(ptr) bind(c)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_associated
   use sim_lw4cuf_class
   implicit none
   type(c_ptr), value :: ptr
   
   !TODO: verify if this can be polymorphic
   type(sim_lw4cuf), pointer :: simp 
   
   if(c_associated(ptr)) then
      call c_f_pointer(ptr, simp)
      print *, "Streaming BW (GB/s): ", ((sizeof(simp%f1) + sizeof(simp%f2))/1.0D9)/(simp%stm/simp%nsteps) 
      print *, "Collision BW (GB/s): ", ((sizeof(simp%f1) + sizeof(simp%f2))/1.0D9)/(simp%ctm/simp%nsteps) 
      DEALLOCATE(simp)
   end if

end subroutine


function c_lw4cuf_norm(nx,ny,u,ua) bind(c)
   use, intrinsic :: iso_c_binding, only: c_int, c_double
   use sim_lw4cuf_class
   integer(c_int), value :: nx, ny
   real(c_double), intent(in) :: u(nx,ny), ua(nx,ny)
   real(c_double) :: c_lw4cuf_norm
   c_lw4cuf_norm = norm2(u - ua) / norm2(ua)
end function


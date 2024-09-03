module sim_class

use kinds, only: dp => wp
implicit none

!> Abstact LBM simulation class; assumes a regular-grid is used.
type, abstract :: sim
   integer :: nx, ny
contains
   procedure(sim_init), deferred :: init
   procedure(sim_step), deferred :: step
   procedure(sim_vars), deferred :: vars
end type

abstract interface
   subroutine sim_init(this,nx,ny,dt,rho,u,sigma)
      import dp, sim
      class(sim), intent(out) :: this
      integer, intent(in) :: nx, ny
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: rho(nx,ny), u(nx,ny,2)
      real(dp), intent(in), optional :: sigma(nx,ny,3)
   end subroutine
   subroutine sim_step(this,omega)
      import dp, sim
      class(sim), intent(inout) :: this
      real(dp), intent(in) :: omega
   end subroutine
   subroutine sim_vars(this,rho,u)
      import dp, sim
      class(sim), intent(in) :: this
      real(dp), intent(out), contiguous :: rho(:,:), u(:,:,:) 
   end subroutine
end interface


!> Abstact LBM simulation class; explicit coordinates
type, abstract :: msim
contains
   procedure(msim_init), deferred :: init
   procedure(msim_step), deferred :: step
   procedure(msim_vars), deferred :: vars
end type

abstract interface
   subroutine msim_init(this,n,xc,yc,dt,rho,u,sigma)
      import dp, msim
      class(msim), intent(out) :: this
      integer, intent(in) :: n
      real(dp), intent(in) :: xc(n), yc(n), dt
      real(dp), intent(in) :: rho(n), u(n,2)
      real(dp), intent(in), optional :: sigma(n,3)
   end subroutine
   subroutine msim_step(this,omega)
      import dp, msim
      class(msim), intent(inout) :: this
      real(dp), intent(in) :: omega
   end subroutine
   subroutine msim_vars(this,rho,u)
      import dp, msim
      class(msim), intent(in) :: this
      real(dp), intent(out), contiguous :: rho(:), u(:,:) 
   end subroutine
end interface

end module

module lbm_primitives

   use sim_class, only: wp => dp
   implicit none
   private

   public :: lbm_grid
   public :: lbm_macros
   public :: lbm_collide_and_stream_fused
   public :: lbm_eqinit
   public :: lbm_neqinit

   public :: lbm_periodic_bc_push
   public :: lbm_periodic_bc_pull

   public :: fluid_state, csqr

   interface lbm_eqinit
      module procedure :: lbm_eqinit_fields
#if !defined(__NVCOMPILER)
      module procedure :: lbm_eqinit_callback
#endif
   end interface

   interface lbm_neqinit
      module procedure :: lbm_neqinit_fields
      module procedure :: lbm_neqinit_callback
   end interface

   type :: fluid_state
      real(wp) :: p, u, v
      real(wp) :: sxx, sxy, syy
   end type

#if defined(__GFORTRAN__) || defined(__NVCOMPILER)
   type :: lbm_grid
      integer :: nx, ny
      real(wp), allocatable :: f1(:,:,:), f2(:,:,:)
      real(wp), allocatable :: rho(:,:), u(:,:), v(:,:)
      real(wp) :: omega
   end type
#else
   type :: lbm_grid(nx,ny)
      integer, len :: nx, ny
      real(wp) :: f1(0:nx+1,0:ny+1,0:8), f2(0:nx+1,0:ny+1,0:8)
      real(wp) :: rho(nx,ny)
      real(wp) :: u(nx,ny), v(nx,ny)
      real(wp) :: omega
   end type
#endif

   real(wp), parameter :: w0 = 4.0_wp / 9.0_wp, &
                          ws = 1.0_wp / 9.0_wp, &
                          wd = 1.0_wp / 36.0_wp

   real(wp), parameter :: w(0:8) = [w0,ws,ws,ws,ws,wd,wd,wd,wd]

   real(wp), parameter :: csqr = 1.0_wp/3.0_wp
   real(wp), parameter :: invcsqr = 1.0_wp/csqr

   real(wp), parameter :: rho0 = 1.0_wp

   ! Double-distribution shifting
   logical, parameter :: ddf_shift = .true.

   integer, parameter :: cx(0:8) = [0,1,0,-1,0,1,-1,-1,1]
   integer, parameter :: cy(0:8) = [0,0,1,0,-1,1,1,-1,-1]

   ! Hermite polynomials
   real(wp), parameter :: hx(0:8) = real(cx,wp)
   real(wp), parameter :: hy(0:8) = real(cy,wp)
   real(wp), parameter :: hxx(0:8) = hx**2 - csqr
   real(wp), parameter :: hxy(0:8) = hx*hy
   real(wp), parameter :: hyy(0:8) = hy**2 - csqr
   real(wp), parameter :: hxxy(0:8) = hxx*hy
   real(wp), parameter :: hyyx(0:8) = hyy*hy
   real(wp), parameter :: hxxyy(0:8) = hxx*hyy

contains

   subroutine lbm_macros(nx,ny,nhalo,fsrc,rho,u,v)
      integer, intent(in) :: nx, ny, nhalo
      real(wp), intent(in) :: fsrc(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8)
      real(wp), intent(out) :: rho(nx,ny), u(nx,ny), v(nx,ny)

      real(wp) :: f(0:8)
      integer :: i, j, k

      !$omp parallel do collapse(2) default(private) shared(fsrc,rho,u,v,nx,ny)
      do j = 1, ny
         do i = 1, nx

            ! Gather PDFs
            do k = 0, 8
               f(k) = fsrc(i,j,k)
            end do

            ! density
            rho(i,j) = (((f(5) + f(7)) + (f(6) + f(8))) + &
                      ((f(1) + f(3)) + (f(2) + f(4)))) + f(0)

            if (ddf_shift) rho(i,j) = rho(i,j) + rho0

            ! velocity
            u(i,j) = (((f(5) - f(7)) + (f(8) - f(6))) + (f(1) - f(3))) / rho(i,j)
            v(i,j) = (((f(5) - f(7)) + (f(6) - f(8))) + (f(2) - f(4))) / rho(i,j)

         end do
      end do
      !$omp end parallel do

   end subroutine

   subroutine lbm_eqinit_fields(nx,ny,nhalo,f,p,u,v)
      integer, intent(in) :: nx, ny, nhalo
      real(wp), intent(out) :: f(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,0:8)
      real(wp), intent(in) :: p(nx,ny), u(nx,ny), v(nx,ny)


      real(wp) :: rho
      integer :: i, j

      !$omp parallel do collapse(2) private(rho) shared(nx,ny,f,p,u,v)
      do j = 1, ny
         do i = 1, nx
            rho = rho0 + p(i,j)/csqr
            f(i,j,:) = equilibrium(rho,u(i,j),v(i,j))
         end do
      end do
      !$omp end parallel do

   end subroutine

   subroutine lbm_neqinit_fields(nx,ny,f,p,u,v,sxx,sxy,syy,taulb)
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: f(0:nx+1,0:ny+1,0:8)
      real(wp), intent(in) :: p(nx,ny), u(nx,ny), v(nx,ny), taulb
      real(wp), intent(in) :: sxx(nx,ny), sxy(nx,ny), syy(nx,ny)


      integer, parameter :: cx(0:8) = [0,1,0,-1,0,1,-1,-1,1]
      integer, parameter :: cy(0:8) = [0,0,1,0,-1,1,1,-1,-1]
      real(wp) :: rho
      integer :: i, j, k
      !call lbm_eqinit_fields(nx,ny,f,p,u,v)

      !$omp parallel do collapse(2) default(shared) private(rho)
      do j = 1, ny
         do i = 1, nx
            rho = rho0 + p(i,j)/csqr
            f(i,j,:) = equilibrium(rho,u(i,j),v(i,j))

            do k = 0, 8
               f(i,j,k) = f(i,j,k) - taulb*rho*3*w(k)*(cx(k)**2 * sxx(i,j) &
                  + cy(k)**2*syy(i,j) + 2*cx(k)*cy(k)*sxy(i,j))
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine

#if defined(__NVCOMPILER)
   #warning "No callback initialization"
#else
   subroutine lbm_eqinit_callback(nx,ny,f,fstate,nu)
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: f(0:nx+1,0:ny+1,0:8)
      real(wp), intent(in) :: nu
      interface
         pure function fstate(x,y)
            import wp, fluid_state
            real(wp), intent(in) :: x, y
            type(fluid_state) :: fstate
         end function
      end interface

      real(wp) :: x, y, rho
      type(fluid_state) :: state
      integer :: i, j

      !$omp parallel do collapse(2) &
      !$omp default(private) shared(nx,ny,f,nu)
      do j = 1, ny
         do i = 1, nx

            x = (i-1) + 0.5_wp
            y = (j-1) + 0.5_wp

            state = fstate(x,y)
            rho = rho0 + state%p/csqr
#if 0
            f(i,j,:) = equilibrium(rho,state%u,state%v)
#else
            block
               integer :: k
               real(wp) :: ux, uy, uxx, uxy, uyy, uxxy, uyyx, uxxyy
               real(wp) :: a1xx, a1xy, a1yy, a1xxy, a1yyx, a1xxyy
               real(wp), parameter :: invcsqr = 1._wp/csqr
               real(wp), parameter :: invcsqr2 = invcsqr**2/2._wp
               real(wp), parameter :: invcsqr3 = invcsqr2*invcsqr/3._wp
               real(wp), parameter :: invcsqr4 = invcsqr3*invcsqr/4._wp
               real(wp) :: feq(0:8), fneq(0:8)

               ux = state%u
               uy = state%v
               uxx = ux*ux
               uxy = ux*uy
               uyy = uy*uy 
               uxxy = uy*uxx
               uyyx = ux*uyy
               uxxyy = uy*uxxy


               do k = 0, 8
                  feq(k) = w(k)*rho*(1.0_wp + invcsqr*(hx(k)*ux + hy(k)*uy) + &
                                              invcsqr2*(hxx(k)*uxx + 2._wp*hxy(k)*uxy + hyy(k)*uyy) + &
                                              3._wp*invcsqr3*(hxxy(k)*uxxy + hyyx(k) * uyyx) + &
                                              6._wp*invcsqr4*(hxxyy(k)*uxxyy))
               end do

               a1xx = 1._wp*rho*nu*state%sxx
               a1xy = 1._wp*rho*nu*state%sxy
               a1yy = 1._wp*rho*nu*state%syy

               a1xxy = 2._wp*ux*a1xy + uy*a1xx
               a1yyx = 2._wp*uy*a1xy + ux*a1yy
               a1xxyy = 2._wp*(ux*a1yyx + uy*a1xxy) - ux**2*a1yy - uy**2*a1xx - 4._wp*ux*uy*a1xy

               do k = 0, 8
                 fneq(k) =  w(k) * (invcsqr2*(hxx(k)*a1xx + 2._wp*hxy(k)*a1xy + hyy(k)*a1yy) + &
                                    3._wp*invcsqr3*(hxxy(k)*a1xxy + hyyx(k) * a1yyx) + &
                                    6._wp*invcsqr4*(hxxyy(k)*a1xxyy))
               end do

               f(i,j,:) = feq + fneq

            end block
#endif
         end do
      end do
      !$omp end parallel do

   end subroutine
#endif

   subroutine lbm_neqinit_callback(nx,ny,f,fstate,taulb)
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: f(0:nx+1,0:ny+1,0:8)
      interface
         pure function fstate(x,y)
            import wp, fluid_state
            real(wp), intent(in) :: x, y
            type(fluid_state) :: fstate
         end function
      end interface
      real(wp), intent(in) :: taulb

      integer, parameter :: cx(0:8) = [0,1,0,-1,0,1,-1,-1,1]
      integer, parameter :: cy(0:8) = [0,0,1,0,-1,1,1,-1,-1]
      real(wp) :: rho, x, y
      integer :: i, j, k
      type(fluid_state) :: state
      !call lbm_eqinit_fields(nx,ny,f,p,u,v)

      !$omp parallel do collapse(2) default(private) shared(nx,ny,f,taulb)
      do j = 1, ny
         do i = 1, nx
            x = (i-1) + 0.5_wp
            y = (j-1) + 0.5_wp

            state = fstate(x,y)
            rho = rho0 + state%p/csqr
            f(i,j,:) = equilibrium(rho,state%u,state%v)

            do k = 0, 8
               f(i,j,k) = f(i,j,k) - taulb*rho*3*w(k)*(cx(k)**2 * state%sxx &
                  + cy(k)**2*state%syy + 2*cx(k)*cy(k)*state%sxy)
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine

   pure function equilibrium(rho,ux,uy) result(feq)
      real(wp), intent(in) :: rho, ux, uy
      real(wp) :: feq(0:8)

      real(wp) :: uxx, uyy,  uxpy, uxmy
      real(wp) :: indp

      uxx = ux*ux
      uyy = uy*uy

      uxpy = ux + uy
      uxmy = ux - uy
   
      if (ddf_shift) then
         !
         ! DDF Shifting (for better accuracy)
         !
         indp = -1.5_wp * (uxx + uyy)
         
         feq(0) = w0*rho*(indp)
         feq(1) = ws*rho*(indp + 3.0_wp*ux + 4.5_wp*uxx)
         feq(2) = ws*rho*(indp + 3.0_wp*uy + 4.5_wp*uyy)
         feq(3) = ws*rho*(indp - 3.0_wp*ux + 4.5_wp*uxx)
         feq(4) = ws*rho*(indp - 3.0_wp*uy + 4.5_wp*uyy)

         feq(5) = wd*rho*(indp + 3.0_wp*uxpy + 4.5_wp*uxpy*uxpy)
         feq(7) = wd*rho*(indp - 3.0_wp*uxpy + 4.5_wp*uxpy*uxpy)
         
         feq(6) = wd*rho*(indp - 3.0_wp*uxmy + 4.5_wp*uxmy*uxmy)
         feq(8) = wd*rho*(indp + 3.0_wp*uxmy + 4.5_wp*uxmy*uxmy)

         feq = feq + w*(rho - rho0)
      else
         indp = 1.0_wp - 1.5_wp * (uxx + uyy)

         feq(0) = w0*rho*(indp)
         feq(1) = ws*rho*(indp + 3.0_wp*ux + 4.5_wp*uxx)
         feq(2) = ws*rho*(indp + 3.0_wp*uy + 4.5_wp*uyy)
         feq(3) = ws*rho*(indp - 3.0_wp*ux + 4.5_wp*uxx)
         feq(4) = ws*rho*(indp - 3.0_wp*uy + 4.5_wp*uyy)

         feq(5) = wd*rho*(indp + 3.0_wp*uxpy + 4.5_wp*uxpy*uxpy)
         feq(7) = wd*rho*(indp - 3.0_wp*uxpy + 4.5_wp*uxpy*uxpy)
         
         feq(6) = wd*rho*(indp - 3.0_wp*uxmy + 4.5_wp*uxmy*uxmy)
         feq(8) = wd*rho*(indp + 3.0_wp*uxmy + 4.5_wp*uxmy*uxmy)
      end if

   end function


   subroutine lbm_collide_and_stream_fused(nx,ny,fsrc,fdst,omega)
      implicit none
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: fsrc(0:nx+1,0:ny+1,0:8)
      real(wp), intent(out) :: fdst(0:nx+1,0:ny+1,0:8)
      real(wp), intent(in) :: omega

      real(wp) :: f(0:8), feq(0:8), rho, ux, uy, irho
      real(wp) :: indp, uxx, uyy, uxpy, uxmy

      integer :: i, j, k

      logical, parameter :: push = .true.

#if defined(__NVCOMPILER)
      !$omp parallel do collapse(2) &
#else
      !$omp parallel do simd collapse(2) schedule(simd: static) &
#endif
      !$omp default(private) shared(nx,ny,fsrc,fdst,omega)
      do j = 1, ny
         do i = 1, nx

            ! Gather PDFs
            if (push) then
               !GCC$ UNROLL 9
               do k = 0, 8
                  f(k) = fsrc(i,j,k)
               end do
            else
               f(0) = fsrc(i  ,j  ,0)
               f(1) = fsrc(i-1,j  ,1)
               f(2) = fsrc(i  ,j-1,2)
               f(3) = fsrc(i+1,j  ,3)
               f(4) = fsrc(i  ,j+1,4)
               f(5) = fsrc(i-1,j-1,5)
               f(6) = fsrc(i+1,j-1,6)
               f(7) = fsrc(i+1,j+1,7)
               f(8) = fsrc(i-1,j+1,8)
            end if

            ! density
            rho = (((f(5) + f(7)) + (f(6) + f(8))) + &
                   ((f(1) + f(3)) + (f(2) + f(4)))) + f(0)
            
            if (ddf_shift) rho = rho + rho0

            irho = 1.0_wp/rho

            ! velocity
            ux = (((f(5) - f(7)) + (f(8) - f(6))) + (f(1) - f(3))) / rho
            uy = (((f(5) - f(7)) + (f(6) - f(8))) + (f(2) - f(4))) / rho

            uxx = ux*ux
            uyy = uy*uy
            
            if (ddf_shift) then
               indp = -1.5_wp * (uxx + uyy)
            else
               indp = 1.0_wp - 1.5_wp * (uxx + uyy)
            end if
            
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

            if (ddf_shift) then
               feq = feq + w*(rho - rho0)
            end if

            ! Scatter PDFs
            if (push) then
               fdst(i  ,j  ,0) = omega*(feq(0) - f(0)) + f(0)
               fdst(i+1,j  ,1) = omega*(feq(1) - f(1)) + f(1)
               fdst(i  ,j+1,2) = omega*(feq(2) - f(2)) + f(2)
               fdst(i-1,j  ,3) = omega*(feq(3) - f(3)) + f(3)
               fdst(i  ,j-1,4) = omega*(feq(4) - f(4)) + f(4)
               fdst(i+1,j+1,5) = omega*(feq(5) - f(5)) + f(5)
               fdst(i-1,j+1,6) = omega*(feq(6) - f(6)) + f(6)
               fdst(i-1,j-1,7) = omega*(feq(7) - f(7)) + f(7)
               fdst(i+1,j-1,8) = omega*(feq(8) - f(8)) + f(8)
            else
               !GCC$ UNROLL 9
               do k = 0, 8
                  fdst(i,j,k) = omega*(feq(k) - f(k)) + f(k)
               end do
            end if

         end do
      end do

   end subroutine

   ! Copy the values from the bulk into the HALO
   subroutine lbm_periodic_bc_pull(nx,ny,fdst)
      implicit none
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: fdst(0:nx+1,0:ny+1,0:8)

      ! 6 2 5
      ! 3 0 1
      ! 7 4 8

   ! Process Edges
 
      ! EAST
      fdst(nx+1, 1:ny, 6) = fdst(1, 1:ny, 6)
      fdst(nx+1, 1:ny, 3) = fdst(1, 1:ny, 3)
      fdst(nx+1, 1:ny, 7) = fdst(1, 1:ny, 7)

      ! WEST
      fdst(0, 1:ny, 5) = fdst(nx, 1:ny, 5)
      fdst(0, 1:ny, 1) = fdst(nx, 1:ny, 1)
      fdst(0, 1:ny, 8) = fdst(nx, 1:ny, 8)

      ! NORTH
      fdst(1:nx, ny+1, 7) = fdst(1:nx, 1, 7)
      fdst(1:nx, ny+1, 4) = fdst(1:nx, 1, 4)
      fdst(1:nx, ny+1, 8) = fdst(1:nx, 1, 8)

      ! SOUTH
      fdst(1:nx, 0, 6) = fdst(1:nx, ny, 6)
      fdst(1:nx, 0, 2) = fdst(1:nx, ny, 2)
      fdst(1:nx, 0, 5) = fdst(1:nx, ny, 5)


   ! Process Corners 
   ! (MUST COME AFTER THE EDGES!)
   !
   ! NW  _________  NE
   !    |         |
   !    | 8     7 |
   !    |         |
   !    | 5     6 |
   !    |_________| 
   ! SW            SE
   !

      ! NORTH-EAST
      fdst(nx+1, ny+1, 7) = fdst(   nx,    ny, 7)

      ! SOUTH-EAST
      fdst(nx+1,  0, 6) = fdst(   ny, 1, 6)

      ! NORTH-WEST
      fdst( 0, ny+1, 8) = fdst(1,   ny, 8)

      ! SOUTH-WEST
      fdst( 0,  0, 5) = fdst(1, 1, 5)

   end subroutine


   !> Copy the values from the HALO layer into the bulk
   subroutine lbm_periodic_bc_push(nx,ny,fdst)
      implicit none
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: fdst(0:nx+1,0:ny+1,0:8)

      ! 6 2 5
      ! 3 0 1
      ! 7 4 8

   ! Process Edges
 
      ! EAST
      fdst(nx, 1:ny, 6) = fdst(0, 1:ny, 6)
      fdst(nx, 1:ny, 3) = fdst(0, 1:ny, 3)
      fdst(nx, 1:ny, 7) = fdst(0, 1:ny, 7)

      ! WEST
      fdst(1, 1:ny, 5) = fdst(nx+1, 1:ny, 5)
      fdst(1, 1:ny, 1) = fdst(nx+1, 1:ny, 1)
      fdst(1, 1:ny, 8) = fdst(nx+1, 1:ny, 8)

      ! NORTH
      fdst(1:nx, ny, 7) = fdst(1:nx, 0, 7)
      fdst(1:nx, ny, 4) = fdst(1:nx, 0, 4)
      fdst(1:nx, ny, 8) = fdst(1:nx, 0, 8)

      ! SOUTH
      fdst(1:nx, 1, 6) = fdst(1:nx, ny+1, 6)
      fdst(1:nx, 1, 2) = fdst(1:nx, ny+1, 2)
      fdst(1:nx, 1, 5) = fdst(1:nx, ny+1, 5)


   ! Process Corners 
   ! (MUST COME AFTER THE EDGES!)
   !
   ! NW  _________  NE
   !    |         |
   !    | 8     7 |
   !    |         |
   !    | 5     6 |
   !    |_________| 
   ! SW            SE
   !

      ! NORTH-EAST
      fdst(nx, ny, 7) = fdst(   0,    0, 7)

      ! SOUTH-EAST
      fdst(nx,  1, 6) = fdst(   0, ny+1, 6)

      ! NORTH-WEST
      fdst( 1, ny, 8) = fdst(nx+1,    0, 8)

      ! SOUTH-WEST
      fdst( 1,  1, 5) = fdst(nx+1, ny+1, 5)

   end subroutine
end module


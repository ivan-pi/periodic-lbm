module lbm_dugks

use sim_class, only: wp => dp

implicit none
private

public :: dugks_step
public :: dugks_stream
public :: dugks_collision
public :: dugks_bc

real(wp), parameter :: cx(0:8) = [0,1,0,-1,0,1,-1,-1,1]
real(wp), parameter :: cy(0:8) = [0,0,1,0,-1,1,1,-1,-1]

real(wp), parameter :: w0 = 4._wp / 9._wp, &
                       ws = 1._wp / 9._wp, &
                       wd = 1._wp / 36._wp
real(wp), parameter :: w(0:8) = [w0,ws,ws,ws,ws,wd,wd,wd,wd]

real(wp), parameter :: rho0 = 1.0_wp

logical, parameter :: use_dugks = .true.

contains

   subroutine dugks_step(nx,ny,fsrc,fdst,tau,dt)
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: fsrc(0:nx+1,0:ny+1,0:8)
      real(wp), intent(inout) :: fdst(0:nx+1,0:ny+1,0:8)
      real(wp), intent(in) :: tau, dt

      real(wp) :: omega

      omega = dt/(tau + dt/2)

      ! Collision for a full time-step
      if (use_dugks) then
         
         ! Create a copy of the field
         fdst(1:nx,1:ny,:) = fsrc(1:ny,1:ny,:)
!         fdst = fsrc
         call dugks_collision(nx,ny,fdst,omega) ! \tilde{f}^{+,n}
         call dugks_collision(nx,ny,fsrc,omega=0.75_wp*omega) ! \bar{f}^{+,n}
!         call dugks_bc(nx,ny,fsrc)

         omega = dt/(4*tau  + dt)
         call dugks_stream(nx,ny,fsrc,fdst,dt,omega,use_dugks)

      else

         ! Bardow's method (kind of ugly however)

         fdst = fsrc
         call dugks_stream(nx,ny,fsrc,fdst,dt,omega,use_dugks)
         call dugks_collision(nx,ny,fdst,omega)
         
      end if

   end subroutine

   pure function equilibrium(rho,ux,uy) result(feq)
      !$omp declare simd simdlen(2)
      !$omp declare simd simdlen(4)
      !$omp declare simd simdlen(8)

      real(wp), intent(in) :: rho, ux, uy
      real(wp) :: feq(0:8)

      real(wp) :: uxx, uyy,  uxpy, uxmy
      real(wp) :: indp

      uxx = ux*ux
      uyy = uy*uy

      ! indp = 1.0_wp - 1.5_wp * (uxx + uyy)

      ! feq(0) = w0*rho*(indp)
      ! feq(1) = ws*rho*(indp + 3.0_wp*ux + 4.5_wp*uxx)
      ! feq(2) = ws*rho*(indp + 3.0_wp*uy + 4.5_wp*uyy)
      ! feq(3) = ws*rho*(indp - 3.0_wp*ux + 4.5_wp*uxx)
      ! feq(4) = ws*rho*(indp - 3.0_wp*uy + 4.5_wp*uyy)

      ! feq(5) = wd*rho*(indp + 3.0_wp*uxpy + 4.5_wp*uxpy*uxpy)
      ! feq(7) = wd*rho*(indp - 3.0_wp*uxpy + 4.5_wp*uxpy*uxpy)
      
      ! feq(6) = wd*rho*(indp - 3.0_wp*uxmy + 4.5_wp*uxmy*uxmy)
      ! feq(8) = wd*rho*(indp + 3.0_wp*uxmy + 4.5_wp*uxmy*uxmy)


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


   subroutine dugks_stream(nx,ny,fsrc,fdst,dt,omega,dugks)
      implicit none
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: fsrc(0:nx+1,0:ny+1,0:8) ! \bar{f}^{+,n}
      real(wp), intent(inout) :: fdst(0:nx+1,0:ny+1,0:8) ! \tilde{f}^{+,n}
      real(wp), intent(in) :: dt, omega
      logical, intent(in) :: dugks

      ! Streaming variables
      real(wp) :: cxk, cyk
      real(wp) :: fc, fn, fs, fw, fe, fne, fnw, fse, fsw

      integer :: i, j, k, ii
      integer,parameter :: bx = 200

      ! PDFs at the cell faces
      real(wp) :: cfw(0:8), cfn(0:8), cfe(0:8), cfs(0:8)
      
      real(wp) :: rw, uxw, uyw
      real(wp) :: re, uxe, uye
      real(wp) :: rn, uxn, uyn
      real(wp) :: rs, uxs, uys

      real(wp) :: fluxW, fluxN, fluxE, fluxS
      
      ! regular derivative parameters
      real(wp), parameter :: p2 = 0.5_wp, &
                             p4 = 0.25_wp, &
                             p8 = 0.125_wp

      real(wp) :: vx(0:8), vy(0:8)

      vx = dt * cx
      vy = dt * cy

!      do ii = 1, nx, bx
      !$omp parallel do collapse(2) default(private) shared(nx,ny,fsrc,fdst,omega,dugks,vx,vy)
      do j = 1, ny
         do i = 1, nx
!         do i = ii, min(ii+bx-1,nx)

            !
            ! Interpolation at the cell faces
            ! (all values are needed due to the macros calculation)
            !
            rw = 0; uxw = 0; uyw = 0;
            re = 0; uxe = 0; uye = 0;
            rn = 0; uxn = 0; uyn = 0;
            rs = 0; uxs = 0; uys = 0;

            do k = 0, 8

               fc  = fsrc(i  , j, k)
               
               fe  = fsrc(i+1, j, k)
               fn  = fsrc(i, j+1, k)
               fw  = fsrc(i-1, j, k)
               fs  = fsrc(i, j-1, k)

               fne = fsrc(i+1, j+1, k)
               fnw = fsrc(i-1, j+1, k)
               fsw = fsrc(i-1, j-1, k)
               fse = fsrc(i+1, j-1, k)

               ! WEST cell face
               cfw(k) = p2*(fc + fw) - p2*vx(k)*(fc - fw) - &
                          p8*vy(k)*(fnw + fn - fsw - fs)

               rw =   rw + cfw(k)
               uxw = uxw + cx(k)*cfw(k)
               uyw = uyw + cy(k)*cfw(k)

               ! NORTH cell face
               cfn(k) = p2*(fc + fn) - p2*vy(k)*(fn - fc) - &
                          p8*vx(k)*(fne + fe - fnw - fw)

               rn =   rn + cfn(k)
               uxn = uxn + cx(k)*cfn(k)
               uyn = uyn + cy(k)*cfn(k)

               ! EAST cell face
               cfe(k) = p2*(fc + fe) - p2*vx(k)*(fe - fc) - &
                          p8*vy(k)*(fne + fn - fse - fs)

               re =   re + cfe(k)
               uxe = uxe + cx(k)*cfe(k)
               uye = uye + cy(k)*cfe(k)

               ! SOUTH cell face
               cfs(k) = p2*(fc + fs) - p2*vy(k)*(fc - fs) - &
                          p8*vx(k)*(fse + fe - fsw - fw)

               rs =   rs + cfs(k)
               uxs = uxs + cx(k)*cfs(k)
               uys = uys + cy(k)*cfs(k)

            end do

            re = re + rho0
            rw = rw + rho0
            rs = rs + rho0
            rn = rn + rho0

            !
            ! Update at the edge 
            ! ( In principle we could save some work here )
            !
            ! Step 4.
            !
            ! At each cell interface, ρ^{n+1/2} and u^{n+1/2} are calculated 
            ! from \tilde{f}^{n+1/2} using Eq. (21)

            ! if (dugks) then
            !    call bgk(cfw,omega)
            !    call bgk(cfe,omega)
            !    call bgk(cfn,omega)
            !    call bgk(cfs,omega)
            ! end if

            cfe = cfe + omega*(equilibrium(re,uxe/re,uye/re) - cfe)
            cfw = cfw + omega*(equilibrium(rw,uxw/rw,uyw/rw) - cfw)
            cfn = cfn + omega*(equilibrium(rn,uxn/rn,uyn/rn) - cfn)
            cfs = cfs + omega*(equilibrium(rs,uxs/rs,uys/rs) - cfs)

            do k = 1, 8
            
               fluxE = cfe(k)
               fluxW = cfw(k)
               fluxN = cfn(k)
               fluxS = cfs(k)

               fdst(i,j,k) = fdst(i,j,k) - vx(k)*(fluxE - fluxW) &
                                         - vy(k)*(fluxN - fluxS)

            end do

!         end do    ! blocking
      end do
   end do
      !$omp end parallel do


   contains

      subroutine bgk(f,omega)
         real(wp), intent(inout) :: f(0:8)
         real(wp), intent(in) :: omega

         real(wp) :: rho, ux, uy, fneq(0:8)
         real(wp), parameter :: rho0 = 1.0_wp

         rho = ((((f(5) + f(7)) + (f(6) + f(8))) + &
                      ((f(1) + f(3)) + (f(2) + f(4)))) + f(0)) + rho0

         ux = (((f(5) - f(7)) + (f(8) - f(6))) + (f(1) - f(3))) / rho
         uy = (((f(5) - f(7)) + (f(6) - f(8))) + (f(2) - f(4))) / rho

         fneq = equilibrium(rho, ux, uy) - f
         f = f + omega*fneq

      end subroutine

   end subroutine


   subroutine dugks_collision(nx,ny,pdf,omega)
      implicit none
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: pdf(0:nx+1,0:ny+1,0:8)
      real(wp), intent(in) :: omega

      ! Collision variables
      real(wp) :: rho, irho, ux, uy, fneq(0:8), f(0:8)
      
      integer :: i, j, k

      !
      ! Collision
      !

      ! Here we perform collision on the HALO too!

      !$omp parallel do collapse(2) private(f,k,rho,irho,ux,uy,fneq)
      do j = 0, ny+1
         do i = 0, nx+1

            ! Gather PDFs
            do k = 0, 8
               f(k) = pdf(i,j,k)
            end do

            ! density
            rho = f(0) + (((f(5) + f(7)) + (f(6) + f(8))) + &
                   ((f(1) + f(3)) + (f(2) + f(4)))) + rho0

            irho = 1.0_wp/rho

            ! velocity
            ux = (((f(5) - f(7)) + (f(8) - f(6))) + (f(1) - f(3))) * irho
            uy = (((f(5) - f(7)) + (f(6) - f(8))) + (f(2) - f(4))) * irho

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


   subroutine dugks_bc(nx,ny,fdst)
      implicit none
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: fdst(0:nx+1,0:ny+1,0:8)

      ! In the Lax-Wendroff method, we need to couple all of the 
      ! PDFs and not just the incoming ones!

      ! SOUTH HALO
      fdst(1:nx,   0,:) = fdst(1:nx,ny,:)

      ! NORTH HALO
      fdst(1:nx,ny+1,:) = fdst(1:nx, 1,:)

      ! WEST HALO
      fdst(   0,1:ny,:) = fdst(nx,1:ny,:)

      ! EAST HALO
      fdst(nx+1,1:ny,:) = fdst( 1,1:ny,:)

      !
      ! CORNERS
      !
      fdst(0,0,:) = fdst(nx, ny, :) 
      fdst(nx+1,ny+1,:) = fdst(1, 1,:) 

      fdst(0, ny+1, :) = fdst(nx,1, :) 
      fdst(nx+1,0, :) = fdst(1,ny, :)


   end subroutine

end module



!> Module implementing standard LBM method
module sim_dugks_class
   
   use sim_class, only: sim, dp
   
   use lbm_primitives, only: &
      lbm_grid, lbm_eqinit, lbm_macros

   use lbm_dugks, only: dugks_step_ => dugks_step, dugks_bc

   implicit none

   type, extends(sim) :: sim_dugks
      type(lbm_grid(:,:)), allocatable :: grid
      integer :: flip = 0
      real(dp) :: dt
   contains
      procedure :: init => dugks_init
      procedure :: step => dugks_step
      procedure :: vars => dugks_vars
   end type

contains

   subroutine dugks_init(this,nx,ny,dt,rho,u,sigma)
      class(sim_dugks), intent(out) :: this
      integer, intent(in) :: nx, ny
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: rho(nx,ny), u(nx,ny,2)    ! rho, (u, v)
      real(dp), intent(in), optional :: sigma(nx,ny,3)  ! (sxx, sxy, syy)

      ! Store the time-step
      this%dt = dt

      ! PDTs for the win!
      allocate(lbm_grid(nx,ny) :: this%grid)
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
         call lbm_eqinit(nx,ny,this%grid%f1,rho,u(:,:,1),u(:,:,2))
      end if

   end subroutine

   subroutine dugks_step(this,omega)
      class(sim_dugks), intent(inout) :: this
      real(dp), intent(in) :: omega

      real(dp) :: tau

      ! TODO: probably this tau is wrong, and we need the full formula with dt
      tau = 1.0_dp/omega

      associate( &
         nx => this%grid%nx, ny => this%grid%ny, &
         f1 => this%grid%f1, f2 => this%grid%f2)

         select case(this%flip)
         case(0)
            call dugks_step_(nx,ny,f1,f2,this%dt,tau)
            call dugks_bc(nx,ny,f2)
         case(1)
            call dugks_step_(nx,ny,f2,f1,this%dt,tau)
            call dugks_bc(nx,ny,f1)
         end select

         this%flip = 1 - this%flip

      end associate

   end subroutine

   subroutine dugks_vars(this,rho,u)
      class(sim_dugks), intent(in) :: this
      real(dp), intent(out), contiguous :: rho(:,:), u(:,:,:)
      associate(grid => this%grid)
         select case(this%flip)
         case(0)
            call lbm_macros(grid%nx,grid%ny,grid%f1,rho,u(:,:,1),u(:,:,2))
         case(1)
            call lbm_macros(grid%nx,grid%ny,grid%f2,rho,u(:,:,1),u(:,:,2))
         end select
      end associate
   end subroutine

end module
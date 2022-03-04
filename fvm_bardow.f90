module precision

   implicit none
   private

   public :: wp

   integer, parameter :: sp = kind(1.0)
   integer, parameter :: dp = kind(1.0d0)

   integer, parameter :: wp = dp

end module

module interp

   use precision, only: wp
   implicit none

contains 


   !   u1 ------- u2
   !   |          |
   !   |          |
   !   |          |
   !   |          |
   !   u3 ---o--- u4
   !   |    /| ry |
   !   |   @--    |
   !   |    rx    |
   !   |          |
   !   u5 ------- u6

   function interp_vbox6(u1,u2,u3,u4,u5,u6,rx,ry) result(uinterp)
      real(wp), intent(in) :: u1, u2, u3, u4, u5, u6
      real(wp), intent(in) :: rx, ry

      real(wp) :: uinterp

      real(wp) :: tx
      real(wp) :: v1, v2, v3
      real(wp) :: h1, h2, h3

      tx = 0.5_wp - rx

      ! interpolate linearly in first direction
      v1 = u1 + (u2 - u1)*tx
      v2 = u3 + (u4 - u3)*tx
      v3 = u5 + (u6 - u5)*tx

      ! interpolate quadratically along second direction
      
      h1 = 0.5_wp*ry*(ry + 1)
      h2 = (ry + 1)*(ry - 1)
      h3 = 0.5_wp*ry*(ry - 1)

      uinterp = h1*v1 + h2*v2 + h3*v3

   end function


   !   u5 ------- u3 ------- u1
   !   |          |          |
   !   |          |          |
   !   |          o          |
   !   |        / |          |
   !   |       @--|          |
   !   u6 ------- u4 ------- u2

   function interp_hbox6(u1,u2,u3,u4,u5,u6,rx,ry) result(uinterp)
      real(wp), intent(in) :: u1, u2, u3, u4, u5, u6
      real(wp), intent(in) :: rx, ry

      real(wp) :: uinterp

      real(wp) :: ty
      real(wp) :: v1, v2, v3
      real(wp) :: h1, h2, h3

      ty = 0.5_wp - ry

      ! interpolate linearly in first direction
      v1 = u2 + (u1 - u2)*ty
      v2 = u4 + (u3 - u4)*ty
      v3 = u6 + (u5 - u6)*ty

      ! interpolate quadratically along second direction
      h1 = 0.5_wp*rx*(rx + 1)
      h2 = (rx + 1)*(rx - 1)
      h3 = 0.5_wp*rx*(rx - 1)

      uinterp = h1*v1 + h2*v2 + h3*v3

   end function

end module

module fvm_bardow

   use precision
   implicit none
   private


   public :: wp
   public :: lattice_grid
   

   public :: alloc_grid, dealloc_grid

   public :: set_properties
   public :: perform_step
   public :: update_macros
   public :: output_grid_txt
   public :: set_pdf_to_equilibrium


   public :: collide_bgk
   public :: stream_fdm_bardow
   public :: stream_fvm_bardow
   public :: stream_fdm_sofonea

   character(*), parameter :: FMT_REAL_SP = '(es15.8e2)'

   type :: lattice_grid

      integer :: nx, ny

      real(wp), allocatable :: f(:,:,:,:)
      real(wp), allocatable :: rho(:,:)
      real(wp), allocatable :: ux(:,:), uy(:,:)
      
      real(wp) :: nu, dt
      real(wp) :: omega
      real(wp) :: csqr
      
      integer :: iold, inew
      
      procedure(collision_interface), pointer, pass(grid) :: collision => null()
      procedure(streaming_interface), pointer, pass(grid) :: streaming => null()

      character(len=:), allocatable :: filename, foldername

      character(len=:), allocatable :: logfile
      procedure(gridlog_interface), pointer, pass(grid) :: logger => null()
      integer :: logunit

   end type

   abstract interface
      subroutine collision_interface(grid)
         import lattice_grid
         class(lattice_grid), intent(inout) :: grid
      end subroutine
      subroutine streaming_interface(grid)
         import lattice_grid
         class(lattice_grid), intent(inout) :: grid
      end subroutine
      subroutine gridlog_interface(grid, step)
         import lattice_grid
         class(lattice_grid), intent(in) :: grid
         integer, intent(in) :: step
      end subroutine
   end interface

   real(wp), parameter :: cx(0:8) = [0, 1, 0, -1, 0, 1, -1, -1, 1]
   real(wp), parameter :: cy(0:8) = [0, 0, 1, 0, -1, 1, 1, -1, -1]

   real(wp), parameter :: w0 = 4._wp/9._wp, &
                          ws = 1._wp/9._wp, &
                          wd = 1._wp/36._wp

   real(wp), parameter :: csqr = 1._wp/3._wp
   real(wp), parameter :: invcsqr = 1._wp/csqr

contains

   pure function equilibrium(rho,ux,uy) result(feq)
      real(wp), intent(in) :: rho, ux, uy
      real(wp) :: feq(0:8)

      real(wp) :: uxx, uyy, uxy, usqr, cu

      uxx = ux*ux
      uyy = uy*uy
      uxy = ux*uy

      usqr = uxx + uyy

      feq(0) = w0*rho*(1.0_wp - 1.5_wp*usqr)
      feq(1) = ws*rho*(1.0_wp + 3.0_wp*ux + 4.5_wp*uxx - 1.5_wp*usqr)
      feq(2) = ws*rho*(1.0_wp + 3.0_wp*uy + 4.5_wp*uyy - 1.5_wp*usqr)
      feq(3) = ws*rho*(1.0_wp - 3.0_wp*ux + 4.5_wp*uxx - 1.5_wp*usqr)
      feq(4) = ws*rho*(1.0_wp - 3.0_wp*uy + 4.5_wp*uyy - 1.5_wp*usqr)

      cu = ux + uy
      feq(5) = wd*rho*(1.0_wp + 3.0_wp*cu + 4.5_wp*cu*cu - 1.5_wp*usqr)
      
      cu = uy - ux
      feq(6) = wd*rho*(1.0_wp + 3.0_wp*cu + 4.5_wp*cu*cu - 1.5_wp*usqr)
      
      cu = -uy - ux
      feq(7) = wd*rho*(1.0_wp + 3.0_wp*cu + 4.5_wp*cu*cu - 1.5_wp*usqr)
      
      cu = ux - uy
      feq(8) = wd*rho*(1.0_wp + 3.0_wp*cu + 4.5_wp*cu*cu - 1.5_wp*usqr)

   end function


   subroutine alloc_grid(grid,nx,ny,log)
      type(lattice_grid), intent(out) :: grid
      integer, intent(in) :: nx, ny
      logical, intent(in), optional :: log

      logical :: log_
      character(len=:), allocatable :: logfile_

      grid%nx = nx
      grid%ny = ny

      allocate(grid%f(ny,nx,0:8,2))
      allocate(grid%rho(ny,nx))
      allocate(grid%ux(ny,nx))
      allocate(grid%uy(ny,nx))

      grid%inew = 1
      grid%iold = 2

      log_ = .true.
      if (present(log)) log_ = log

      if (log_) then
         logfile_ = "lattice_grid_log.txt"
         if (allocated(grid%logfile)) then
            logfile_ = grid%logfile
         end if
         open(newunit=grid%logunit,file=logfile_,status='unknown')
      end if

   end subroutine

   subroutine dealloc_grid(grid)
      type(lattice_grid), intent(inout) :: grid
      
      ! close log file
      close(grid%logunit)

      call dealloc_grid_(grid)
   contains
      !> Deallocate all allocatable objects by applying intent(out)
      subroutine dealloc_grid_(grid)
         type(lattice_grid), intent(out) :: grid
         ! Trick to prevent spurious warning
         associate(nx => grid%nx)
         end associate
      end subroutine
   end subroutine

   subroutine set_properties(grid, nu, dt)
      type(lattice_grid), intent(inout) :: grid
      real(wp), intent(in) :: nu, dt

      real(wp) :: tau

      grid%nu = nu
      grid%dt = dt

      grid%csqr = csqr
      tau = invcsqr*nu

      grid%omega = dt/(tau + 0.5_wp*dt)

   end subroutine


   subroutine set_pdf_to_equilibrium(grid)
      type(lattice_grid), intent(inout) :: grid

      real(wp) :: rho_, ux_, uy_
      integer :: x, y

      associate( nx => grid%nx, &
                 ny => grid%ny, &
                 pdf => grid%f(:,:,:,grid%iold), &
                 rho => grid%rho, &
                 ux => grid%ux, &
                 uy => grid%uy)


      do y = 1, ny
         do x = 1, nx

            rho_ = rho(y,x)
            ux_ = ux(y,x)
            uy_ = uy(y,x)

            pdf(y,x,:) = equilibrium(rho_, ux_, uy_)

         end do
      end do

      end associate

   end subroutine

   subroutine perform_step(grid)
      type(lattice_grid), intent(inout) :: grid

      call grid%streaming()
      call grid%collision()

      swap: block
         integer :: itmp
         itmp = grid%iold
         grid%iold = grid%inew
         grid%inew = itmp         
      end block swap

   end subroutine

   subroutine update_macros(grid)
      type(lattice_grid), intent(inout) :: grid

      call update_macros_kernel(grid%nx, grid%nx, &
         grid%f(:,:,:,grid%inew), &
         grid%rho, grid%ux, grid%uy)

   contains

      subroutine update_macros_kernel(nx,ny,f,grho,gux,guy)
         integer, intent(in) :: nx, ny
         real(wp), intent(in) :: f(ny, nx, 0:8)
         real(wp), intent(out), dimension(ny,nx) :: grho, gux, guy

         real(wp) :: rho, invrho, ux, uy, fs(0:8)
         integer :: x, y

         !$omp parallel do collapse(2) default(private) shared(nx,ny,f,rho,gux,guy)
         do x = 1, nx
            do y = 1, ny

               fs = f(y,x,:)

               ! density
               rho = (((fs(5) + fs(7)) + (fs(6) + fs(8))) + &
                      ((fs(1) + fs(3)) + (fs(2) + fs(4)))) + fs(0)

               invrho = 1.0_wp/rho

               ! velocity
               ux = invrho * (((fs(5) - fs(7)) + (fs(8) - fs(6))) + (fs(1) - fs(3)))
               uy = invrho * (((fs(5) - fs(7)) + (fs(6) - fs(8))) + (fs(2) - fs(4)))

               grho(y,x) = rho
               gux(y,x) = ux
               guy(y,x) = uy

            end do
         end do
         !$omp end parallel do

      end subroutine

   end subroutine



   subroutine collide_bgk(grid)
      class(lattice_grid), intent(inout) :: grid

      call bgk_kernel(grid%nx, grid%ny,&
         grid%f(:,:,:,grid%inew), &
         grid%omega)

   contains

      subroutine bgk_kernel(nx,ny,f1,omega)
         integer, intent(in) :: nx, ny
         real(wp), intent(inout) :: f1(ny,nx,0:8)
         real(wp), intent(in) :: omega

         real(wp) :: rho, invrho, ux, uy
         real(wp) :: omegabar
         real(wp) :: fs(0:8), feq(0:8)

         integer :: x, y

         omegabar = 1.0_wp - omega

         !$omp parallel do collapse(2) default(private) shared(f1,omega,omegabar)
         do y = 1, ny
            do x = 1, nx

            ! pull pdfs travelling in different directions
            fs = f1(y,x,:)

            ! density
            rho = (((fs(5) + fs(7)) + (fs(6) + fs(8))) + &
                   ((fs(1) + fs(3)) + (fs(2) + fs(4)))) + fs(0)

            invrho = 1.0_wp/rho

            ! velocity
            ux = invrho * (((fs(5) - fs(7)) + (fs(8) - fs(6))) + (fs(1) - fs(3)))
            uy = invrho * (((fs(5) - fs(7)) + (fs(6) - fs(8))) + (fs(2) - fs(4)))

            ! get equilibrium pdfs
            feq = equilibrium(rho, ux, uy)

            ! collision
            fs = omegabar*fs + omega*feq

            ! push pdfs to destination array
            f1(y,x,:) = fs

            end do
         end do
         !$omp end parallel do

      end subroutine

   end subroutine

   subroutine stream_fvm_bardow(grid)

      use interp, only: interp_hbox6, interp_vbox6

      class(lattice_grid), intent(inout) :: grid

      call fvm_bardow_kernel( &
         grid%nx, grid%ny, &
         grid%f(:,:,:,grid%iold), &
         grid%f(:,:,:,grid%inew), &
         grid%dt)

   contains

      subroutine fvm_bardow_kernel(nx,ny,fold,fnew,dt)
         integer, intent(in) :: nx, ny
         real(wp), intent(in) :: fold(ny,nx,0:8)
         real(wp), intent(inout) :: fnew(ny,nx,0:8)
         real(wp), intent(in) :: dt

         integer :: q, x, y
         integer :: xp1, yp1, xm1, ym1

         real(wp) :: cxq, cyq
         real(wp) :: fc, fe, fn, fw, fs, fne, fnw, fsw, fse
         real(wp) :: cfe, cfn, cfw, cfs

         real(wp), parameter :: p2 = 0.5_wp, p8 = 0.125_wp, p6 = 1._wp/6._wp

         ! copy rest populations
         fnew(:,:,0) = fold(:,:,0)

         do q = 1, 8
         
            cxq = dt*cx(q)
            cyq = dt*cy(q)

         do x = 1, ny

            ! locate neighbor nodes
            xp1 = mod(x, nx) + 1
            xm1 = mod(nx + x - 2, nx) + 1

            do y = 1, nx

               yp1 = mod(y, ny) + 1
               ym1 = mod(ny + y - 2, ny) + 1

               ! read pdf values
               fc  = fold(y  , x  , q)
               fe  = fold(y  , xp1, q)
               fn  = fold(yp1, x  , q)
               fw  = fold(y  , xm1, q)
               fs  = fold(ym1, x  , q)
               fne = fold(yp1, xp1, q)
               fnw = fold(yp1, xm1, q)
               fsw = fold(ym1, xm1, q)
               fse = fold(ym1, xp1, q)

               ! interpolate f at cell faces using
               ! bilinear interpolation
               
               cfw = p2*(fc + fw) - p2*cxq*(fc - fw) - &
                     p8*cyq*(fnw + fn - fsw - fs)

               cfn = p2*(fc + fn) - p2*cyq*(fn - fc) - &
                     p8*cxq*(fne + fe - fnw - fw)

               cfe = p2*(fc + fe) - p2*cxq*(fe - fc) - &
                     p8*cyq*(fne + fn - fse - fs)

               cfs = p2*(fc + fs) - p2*cyq*(fc - fs) - &
                     p8*cxq*(fse + fe - fsw - fw)


               ! Least-squares reconstruction
               !
               !cfw = p2*(fc + fw) - p6*cxq*(fc - fw + fn - fnw + fs - fsw) - &
               !      p8*cyq*(fnw + fn - fsw - fs)

               !cfn = p2*(fc + fn) - p6*cyq*(fn - fc + fne - fe + fnw - fw) - &
               !      p8*cxq*(fne + fe - fnw - fw)

               !cfe = p2*(fc + fe) - p6*cxq*(fe - fc + fne - fn + fse - fs) - &
               !      p8*cyq*(fne + fn - fse - fs)

               !cfs = p2*(fc + fs) - p6*cyq*(fc - fs + fe - fse + fw - fsw) - &
               !      p8*cxq*(fse + fe - fsw - fw)

               ! Lagrangian-tensor interpolation, doesn't work very well however
               !cfw = interp_vbox6(fnw, fn, fw, fc, fsw, fs, -p2*cxq, -p2*cyq)
               !cfe = interp_vbox6(fn, fne, fc, fe, fs, fse, -p2*cxq, -p2*cyq)
               !cfn = interp_hbox6(fne, fe, fn, fc, fnw, fw, -p2*cxq, -p2*cyq)
               !cfs = interp_hbox6(fe, fse, fc, fs, fw, fsw, -p2*cxq, -p2*cyq)

               ! calculate updated value from surface fluxes
               fnew(y, x, q) = fc - cxq*(cfe - cfw) - cyq*(cfn - cfs)

            end do
         end do

         end do

      end subroutine

   end subroutine

   subroutine stream_fdm_bardow(grid)
      class(lattice_grid), intent(inout) :: grid

      call fdm_bardow_kernel( &
         grid%nx, grid%ny, &
         grid%f(:,:,:,grid%iold), &
         grid%f(:,:,:,grid%inew), &
         grid%dt)

   contains

      subroutine fdm_bardow_kernel(nx, ny,fold,fnew,dt)
         integer, intent(in) :: nx, ny
         real(wp), intent(in) :: fold(ny,nx,0:8)
         real(wp), intent(inout) :: fnew(ny,nx,0:8)
         real(wp), intent(in) :: dt

         integer :: x, y, q
         integer :: xp1, yp1, xm1, ym1

         real(wp) :: fc, fe, fn, fw, fs, fne, fnw, fsw, fse

         real(wp) :: dfx, dfy, dfxx, dfxy, dfyy
         real(wp) :: cxq, cyq, cxxq, cyyq, cxyq

         ! regular derivative parameters
         real(wp), parameter :: p2 = 0.5_wp, p4 = 0.25_wp
         
         ! isotropic derivative parameters
         real(wp), parameter :: two_thirds = 2._wp/3._wp, one_sixth = 1._wp/6._wp
         real(wp), parameter :: five_sixths = 10._wp/12._wp, one_twelth = 1._wp/12._wp

         ! copy rest populations
         !$omp parallel default(private) shared(nx,ny,fold,fnew,dt)

         !$omp workshare
         fnew(:,:,0) = fold(:,:,0)
         !$omp end workshare

         do q = 1, 8
         
            cxq = dt*cx(q)
            cyq = dt*cy(q)

            cxxq = cxq*cxq
            cyyq = cyq*cyq
            cxyq = 2*cxq*cyq

            !$omp do
            do x = 1, nx

               ! locate neighbor nodes
               xp1 = mod(x, nx) + 1
               xm1 = mod(nx + x - 2, nx) + 1

               do y = 1, ny

                  yp1 = mod(y, ny) + 1
                  ym1 = mod(ny + y - 2, ny) + 1

                  ! read old pdf values
                  fc  = fold(y  , x  , q)
                  fe  = fold(y  , xp1, q)
                  fn  = fold(yp1, x  , q)
                  fw  = fold(y  , xm1, q)
                  fs  = fold(ym1, x  , q)
                  fne = fold(yp1, xp1, q)
                  fnw = fold(yp1, xm1, q)
                  fsw = fold(ym1, xm1, q)
                  fse = fold(ym1, xp1, q)

                  ! first derivatives
                  dfx = p2*(fe - fw)
                  dfy = p2*(fn - fs)

                  ! second derivatives
                  dfxx = fe - 2*fc + fw
                  dfyy = fn - 2*fc + fs
                  dfxy = p4*(fne - fse - fnw + fsw)

                  !> isotropic derivatives, taken from
                  ! 
                  !   Kumar, A. (2004). Isotropic finite-differences. Journal of
                  !   Computational Physics, 201(1), 109-118, doi:10.1016/j.jcp.2004.05.005
                  !
                  !   the mixed derivative is already isotropic    
                  !
                  !dfx = p2*(one_sixth*(fne - fnw) + two_thirds*(fe - fw) + one_sixth*(fse - fsw))
                  !dfy = p2*(one_sixth*(fne - fse) + two_thirds*(fn - fs) + one_sixth*(fnw - fsw))

                  !dfxx = one_twelth*(fne - 2*fn + fnw) + five_sixths*dfxx + one_twelth*(fse - 2*fs + fsw)
                  !dfyy = one_twelth*(fne - 2*fe + fse) + five_sixths*dfyy + one_twelth*(fnw - 2*fw + fsw)

                  !
                  ! perform streaming
                  !
                  fnew(y,x,q) = fc - cxq*dfx - cyq*dfy + &
                     p2*(cxxq*dfxx + cxyq*dfxy + cyyq*dfyy)

               end do
            end do
            !$omp end do

         end do

         !$omp end parallel
      end subroutine

   end subroutine


   subroutine stream_fdm_sofonea(grid)
      class(lattice_grid), intent(inout) :: grid

      call fdm_sofonea_kernel( &
         grid%nx, grid%ny, &
         grid%f(:,:,:,grid%iold), &
         grid%f(:,:,:,grid%inew), &
         grid%dt)

   contains

      subroutine fdm_sofonea_kernel(nx, ny,fold,fnew,dt)
         integer, intent(in) :: nx, ny
         real(wp), intent(in) :: fold(ny,nx,0:8)
         real(wp), intent(inout) :: fnew(ny,nx,0:8)
         real(wp), intent(in) :: dt

         integer :: x, y
         integer :: xp1, yp1, xm1, ym1

         real(wp) :: du1, du2, fu, fc, fd

         real(wp) :: p2 = 0.5_wp/sqrt(2._wp)

         ! copy rest populations
         fnew(:,:,0) = fold(:,:,0)

         ! Direction 1
         !
         do x = 1, ny
            xp1 = mod(x, nx) + 1
            xm1 = mod(nx + x - 2, nx) + 1
            do y = 1, nx

               fc = fold(y,x  ,1)
               fu = fold(y,xp1,1)
               fd = fold(y,xm1,1)

               du1 = 0.5_wp*(fu - fd)
               du2 = fu - 2*fc + fd
               fnew(y,x,1) = fc + dt*(0.5_wp*dt*du2 - du1)
            end do
         end do

         ! Direction 2
         !
         do x = 1, ny
            do y = 1, nx
               yp1 = mod(y, ny) + 1
               ym1 = mod(ny + y - 2, ny) + 1

               fc = fold(y,  x,2)
               fu = fold(yp1,x,2)
               fd = fold(ym1,x,2)

               du1 = 0.5_wp*(fu - fd)
               du2 = fu - 2*fc + fd
               fnew(y,x,2) = fc + dt*(0.5_wp*dt*du2 - du1)
            end do
         end do

         ! Direction 3
         !
         do x = 1, ny
            xp1 = mod(x, nx) + 1
            xm1 = mod(nx + x - 2, nx) + 1
            do y = 1, nx

               fc = fold(y,x  ,3)
               fu = fold(y,xm1,3)
               fd = fold(y,xp1,3)

               du1 = 0.5_wp*(fu - fd)
               du2 = fu - 2*fc + fd
               fnew(y,x,3) = fc + dt*(0.5_wp*dt*du2 - du1)
            end do
         end do

         ! Direction 4
         !
         do x = 1, ny
            do y = 1, nx
               yp1 = mod(y, ny) + 1
               ym1 = mod(ny + y - 2, ny) + 1

               fc = fold(y,  x,4)
               fu = fold(ym1,x,4)
               fd = fold(yp1,x,4)

               du1 = 0.5_wp*(fu - fd)
               du2 = fu - 2*fc + fd
               fnew(y,x,4) = fc + dt*(0.5_wp*dt*du2 - du1)
            end do
         end do

         ! Direction 5
         !
         do x = 1, ny
            ! locate neighbor nodes
            xp1 = mod(x, nx) + 1
            xm1 = mod(nx + x - 2, nx) + 1
            do y = 1, nx
               yp1 = mod(y, ny) + 1
               ym1 = mod(ny + y - 2, ny) + 1

               fc = fold(y,x,5)
               fu = fold(yp1,xp1,5)
               fd = fold(ym1,xm1,5)

               du1 = p2*(fu - fd)
               du2 = 0.5_wp*(fu - 2*fc + fd)
               fnew(y,x,5) = fc + dt*(0.5_wp*dt*du2 - du1)
            end do
         end do

         ! Direction 6
         !
         do x = 1, ny
            ! locate neighbor nodes
            xp1 = mod(x, nx) + 1
            xm1 = mod(nx + x - 2, nx) + 1
            do y = 1, nx
               yp1 = mod(y, ny) + 1
               ym1 = mod(ny + y - 2, ny) + 1

               fc = fold(y,x,6)
               fu = fold(yp1,xm1,6)
               fd = fold(ym1,xp1,6)

               du1 = p2*(fu - fd)
               du2 = 0.5_wp*(fu - 2*fc + fd)
               fnew(y,x,6) = fc + dt*(0.5_wp*dt*du2 - du1)
            end do
         end do

         ! Direction 7
         !
         do x = 1, ny
            ! locate neighbor nodes
            xp1 = mod(x, nx) + 1
            xm1 = mod(nx + x - 2, nx) + 1
            do y = 1, nx
               yp1 = mod(y, ny) + 1
               ym1 = mod(ny + y - 2, ny) + 1

               fc = fold(y,x,7)
               fu = fold(ym1,xm1,7)
               fd = fold(yp1,xp1,7)

               du1 = p2*(fu - fd)
               du2 = 0.5_wp*(fu - 2*fc + fd)
               fnew(y,x,7) = fc + dt*(0.5_wp*dt*du2 - du1)
            end do
         end do

         ! Direction 8
         !
         do x = 1, ny
            ! locate neighbor nodes
            xp1 = mod(x, nx) + 1
            xm1 = mod(nx + x - 2, nx) + 1
            do y = 1, nx
               yp1 = mod(y, ny) + 1
               ym1 = mod(ny + y - 2, ny) + 1

               fc = fold(y,x,8)
               fu = fold(ym1,xp1,8)
               fd = fold(yp1,xm1,8)

               du1 = p2*(fu - fd)
               du2 = 0.5_wp*(fu - 2*fc + fd)
               fnew(y,x,8) = fc + dt*(0.5_wp*dt*du2 - du1)
            end do
         end do
      end subroutine

   end subroutine

   subroutine output_grid_txt(grid, step)
      type(lattice_grid), intent(in) :: grid
      integer, intent(in), optional :: step

      character(len=64) :: istr
      character(len=:), allocatable :: fullname
      integer :: unit, istat, x, y
      real(wp) :: xx, yy

      istr = ''
      if (present(step)) then
         write(istr,'(I0.9)') step
      end if

      fullname = ''
      if (allocated(grid%foldername)) then
         call execute_command_line("mkdir -p "//trim(grid%foldername), &
            exitstat=istat, wait=.true.)

         if (istat /= 0) then
            write(*,'(A)') "[output_grid] error making directory "//grid%foldername
            error stop
         end if
         fullname = trim(grid%foldername)
      end if

      fullname = fullname//'/'//grid%filename//trim(istr)//'.txt'
      
      open(newunit=unit,file=fullname)

      do x = 1, grid%nx
         xx = (x - 1) + 0.5_wp
         do y = 1, grid%ny
            yy = (y - 1) + 0.5_wp
            write(unit,'(*(1X,'//FMT_REAL_SP//'))') xx, yy, grid%rho(y,x), grid%ux(y,x), grid%uy(y,x)
         end do
         write(unit,*)
      end do

      close(unit)

   end subroutine


end module

module taylor_green

   use precision, only: wp
   use fvm_bardow, only: lattice_grid
   
   implicit none
   private

   public :: taylor_green_t
   public :: pi

   type :: taylor_green_t
      integer :: nx, ny
      real(wp) :: kx, ky
      real(wp) :: umax, nu
      real(wp) :: td
   contains
      procedure :: eval => taylor_green_eval
      procedure :: decay_time => taylor_green_decay_time
   end type

   real(wp), parameter :: pi = 4._wp*atan(1._wp)


   interface taylor_green_t
      module procedure :: taylor_green_t_constructor
   end interface

contains

   function taylor_green_t_constructor(nx,ny,kx,ky,umax,nu) result(this)
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: kx, ky, umax, nu
      type(taylor_green_t) :: this

      this%nx = nx
      this%ny = ny
      this%kx = kx
      this%ky = ky
      this%umax = umax
      this%nu = nu
      this%td = 1._wp/(nu*(kx**2 + ky**2))
   end function

   function taylor_green_decay_time(self) result(tc)
      class(taylor_green_t), intent(in) :: self
      real(wp) :: tc
      tc = self%td
   end function

   subroutine taylor_green_eval(self,t,p,ux,uy,S)
      class(taylor_green_t), intent(in) :: self
      real(wp), intent(in) :: t
      real(wp), intent(out) :: p(:,:), ux(:,:), uy(:,:)
      real(wp), intent(out), optional :: S(self%ny,self%nx,3)

      integer :: x, y
      real(wp) :: xx, yy
      real(wp), parameter :: rho0 = 1.0_wp

      associate(umax=>self%umax, &
                kx => self%kx, &
                ky=>self%ky, &
                td => self%td)

      do x = 1, self%nx
         xx = (x - 1) + 0.5_wp
         do y = 1, self%nx
            yy = (y - 1) + 0.5_wp

            ux(y,x) = -umax*sqrt(ky/kx)*cos(kx*xx)*sin(ky*yy)*exp(-t/td)
            uy(y,x) =  umax*sqrt(kx/ky)*sin(kx*xx)*cos(ky*yy)*exp(-t/td)
            p(y,x) = -0.25_wp*(umax**2)*((ky/kx)*cos(2._wp*kx*xx) + (kx/ky)*cos(2._wp*ky*yy))*exp(-2._wp*t/td)

            if (present(S)) then
              S(y,x,1) = umax*sqrt(kx*ky)*sin(kx*xx)*sin(ky*yy)*exp(-t/td)
              S(y,x,2) = 0.5_wp*umax*(sqrt(kx**3/ky) - sqrt(ky**3/kx))*cos(kx*xx)*cos(ky*yy)*exp(-t/td)
              S(y,x,3) = -S(y,x,1)
            end if
         end do
      end do

      end associate
   end subroutine

!    kykx = ky/kx
!    kxky = kx/ky

!    pfx = -u0*sqrt(kykx)
!    pfy =  u0*sqrt(kxky)

!    do i = 1, n
!      ux(i) = pfx*cos(kx*x(i))*sin(ky*y(i))*exp(-t/tc)
!      uy(i) = pfy*sin(kx*x(i))*cos(ky*y(i))*exp(-t/tc)
!      p(i) = -0.25_wp*u0**2*(kykx*cos(2*kx*x(i)) + kxky*cos(2*ky*y(i)))*exp(-2.0_wp*t/tc)
!    end do

end module

program main

   use fvm_bardow
   use taylor_green, only: taylor_green_t, pi

   implicit none

   integer, parameter :: nx = 64, ny = 64
   integer, parameter :: nprint = 10000
   integer :: step, nsteps

   type(taylor_green_t) :: tg

   type(lattice_grid) :: grid
   real(wp) :: cfl, dt, nu, tau
   real(wp) :: kx, ky, umax, nrm

   real(wp) :: t, tmax

   call alloc_grid(grid, nx, ny)

   grid%filename = "results"
   grid%foldername = "taylor_green"

   grid%collision => collide_bgk
   grid%streaming => stream_fdm_bardow

   grid%logger => my_logger

   ! umax = Mach * cs
   umax = 0.01_wp / sqrt(3._wp)

   ! nu = (umax * L) / Re
   nu = (umax * real(nx,wp)) / 100._wp
   
   ! tau = nu * cs**2
   tau = nu / 3._wp
   dt = 40._wp*tau

   !cfl = 0.1_wp
   !dt = cfl/sqrt(2.0_wp)
   cfl = sqrt(2._wp)*dt
   print *, "dt/tau = ", dt/tau
   print *, "cfl = ", cfl

   call set_properties(grid, nu, dt)

   print *, "omega = ", grid%omega

   ! ---- prepare flow case ----

   kx = 2*pi/real(nx,wp)
   ky = 2*pi/real(ny,wp)

   tg = taylor_green_t(nx,ny,kx,ky,umax,nu)
   print *, "umax = ", umax
   print *, "tc   = ", tg%td
   call write_gnuplot_include()

   tmax = log(2._wp)*tg%decay_time()
   nsteps = int(1.1_wp*tmax/dt)
   print*, "nsteps = ", nsteps

   t = 0._wp
   call apply_initial_condition(tg, grid)

   call output_grid_txt(grid,step=0)
   call grid%logger(step=0)

   time: do step = 1, nsteps

      call perform_step(grid)
      t = t + dt

      if (mod(step,nprint) == 0) then
         print *, step
         call update_macros(grid)
         call output_grid_txt(grid,step)
         call grid%logger(step)
      end if

      if (t >= tmax) then
         call update_macros(grid)
         call output_grid_txt(grid,step)
         call grid%logger(step)
         exit time
      end if
   end do time


   ! calculate average L2-norm
   nrm = calc_L2_norm(tg, grid, t)
   print *, "L2-norm = ", nrm

   call dealloc_grid(grid)

contains

   subroutine apply_initial_condition(case, grid)
      type(taylor_green_t), intent(in) :: case
      type(lattice_grid), intent(inout) :: grid

      real(wp), parameter :: rho0 = 1.0_wp

      call case%eval(t=0.0_wp, &
                   p=grid%rho, &
                   ux=grid%ux, &
                   uy=grid%uy)

      ! convert pressure to lattice density
      grid%rho = grid%rho/grid%csqr + rho0

      call set_pdf_to_equilibrium(grid)

   end subroutine

   subroutine my_logger(grid,step)
      class(lattice_grid), intent(in) :: grid
      integer, intent(in) :: step

      write(grid%logunit, *) step, step*dt, maxval(hypot(grid%ux,grid%uy))
      flush(grid%logunit)

   end subroutine

   subroutine write_gnuplot_include()

      integer :: unit

      open(newunit=unit,file="lattice_grid_log.incl",status='unknown')

      write(unit,*) "dt = ", dt
      write(unit,*) "umax = ", umax
      write(unit,*) "tc = ", tg%td

      close(unit)
   end subroutine


   function calc_L2_norm(case, grid, t) result(nrm)
      type(taylor_green_t), intent(in) :: case
      type(lattice_grid), intent(in) :: grid
      real(wp), intent(in) :: t
      real(wp) :: nrm

      real(wp), allocatable :: pa(:,:), uxa(:,:), uya(:,:)
      real(wp) :: above, below

      allocate(pa, mold=grid%rho) ! not needed
      allocate(uxa, mold=grid%ux)
      allocate(uya, mold=grid%uy)

      call case%eval(t=t, &
                   p=pa, &
                   ux=uxa, &
                   uy=uya)

      above = norm2(hypot(grid%ux-uxa, grid%uy-uya))
      below = norm2(hypot(uxa, uya))
      nrm = above/below

      ! # Code used in the Python version
      ! def L2_error(u,v,ua,va):
      !     return np.sqrt(np.sum((u-ua)**2 + (v-va)**2)/np.sum(ua**2 + va**2))

      !above = sum((grid%ux - uxa)**2 + (grid%uy - uya)**2)
      !below = sum(uxa**2 + uya**2)
      !nrm = sqrt(above/below)

   end function

end program
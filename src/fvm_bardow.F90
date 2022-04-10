module fvm_bardow

   use precision
   use vtk, only: output_vtk_grid_ascii, output_vtk_structuredPoints
   use output_gnuplot, only: output_gnuplot_grid
   
   implicit none
   private


   public :: wp
   public :: lattice_grid
   

   public :: alloc_grid, dealloc_grid

   public :: set_properties
   public :: perform_step, perform_triple_step
   
   public :: update_macros
   
   public :: output_grid_txt, output_vtk

   public :: set_pdf_to_equilibrium

   public :: collide_bgk
   public :: stream_fdm_bardow
   public :: stream_fvm_bardow
   public :: stream_fdm_sofonea

   public :: cx, cy, csqr

   character(*), parameter :: FMT_REAL_SP = '(es15.8e2)'

   type :: lattice_grid

      integer :: nx, ny

      real(wp), allocatable :: f(:,:,:,:)
      !dir$ attributes align: 64:: f

      real(wp), allocatable :: rho(:,:)
      real(wp), allocatable :: ux(:,:), uy(:,:)
      !dir$ attributes align: 64:: rho, ux, uy
      
      real(wp) :: nu, dt, tau
      real(wp) :: omega, trt_magic
      real(wp) :: csqr
      
      integer :: iold, inew, imid
      
      procedure(collision_interface), pointer, pass(grid) :: collision => null()
      procedure(streaming_interface), pointer, pass(grid) :: streaming => null()

      character(len=:), allocatable :: filename, foldername

      character(len=:), allocatable :: logfile
      procedure(gridlog_interface), pointer, pass(grid) :: logger => null()
      integer :: logunit
   contains
      procedure :: set_output_folder
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

   real(wp), parameter :: cx(0:8) = [real(wp) :: 0, 1, 0, -1, 0, 1, -1, -1, 1]
   real(wp), parameter :: cy(0:8) = [real(wp) :: 0, 0, 1, 0, -1, 1, 1, -1, -1]

   real(wp), parameter :: w0 = 4._wp / 9._wp, &
                          ws = 1._wp / 9._wp, &
                          wd = 1._wp / 36._wp

   real(wp), parameter :: csqr = 1._wp/3._wp
   real(wp), parameter :: invcsqr = 1._wp/csqr

contains

   pure function equilibrium(rho,ux,uy) result(feq)
      real(wp), intent(in) :: rho, ux, uy
      real(wp) :: feq(0:8)

      real(wp) :: uxx, uyy, uxy, uxpy, uxmy
      real(wp) :: indp

      uxx = ux*ux
      uyy = uy*uy
      uxy = ux*uy

      indp = 1.0_wp - 1.5_wp * (uxx + uyy)

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

   end function


   subroutine alloc_grid(grid,nx,ny,nf,log)
      type(lattice_grid), intent(out) :: grid
      
      integer, intent(in) :: nx, ny
      integer, intent(in), optional :: nf
      logical, intent(in), optional :: log

      integer :: nf_
      logical :: log_
      character(len=:), allocatable :: logfile_
      integer :: ny_, rmd

      grid%nx = nx
      grid%ny = ny

      ! Round to closest dimension of 16
      ny_ = ny
      rmd = mod(ny_, 16)      
      if (rmd > 0) ny_ = ny_ + (16 - rmd)

      ! Number of pdf fields
      nf_ = 2
      if (present(nf)) nf_ = nf

      allocate(grid%f(ny_,nx,0:8,nf_))
      allocate(grid%rho(ny,nx))
      allocate(grid%ux(ny,nx))
      allocate(grid%uy(ny,nx))

      grid%inew = 1
      grid%iold = 2

      if (nf_ > 2) then
         grid%imid = 3
      else
         grid%imid = -1
      end if

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

   subroutine set_properties(grid, nu, dt, magic)
      type(lattice_grid), intent(inout) :: grid
      real(wp), intent(in) :: nu, dt
      real(wp), optional :: magic

      real(wp) :: tau

      grid%nu = nu
      grid%dt = dt

      grid%csqr = csqr
      tau = invcsqr*nu

      grid%tau = tau
      grid%omega = dt/(tau + 0.5_wp*dt)

      if (present(magic)) then
         grid%trt_magic = magic
      else
         ! bgk by default
         !grid%trt_magic = (2.0_wp - grid%omega)**2 / (4.0_wp * (grid%omega)**2)
         grid%trt_magic = (tau/dt)**2
      end if

      print *, "trt magic = ", grid%trt_magic

   end subroutine


   subroutine set_pdf_to_equilibrium(grid)
      type(lattice_grid), intent(inout) :: grid

      real(wp) :: rho_, ux_, uy_
      integer :: x, y

      associate( nx => grid%nx, &
                 ny => grid%ny, &
                 rho => grid%rho, &
                 ux => grid%ux, &
                 uy => grid%uy, &
                 pdf => grid%f(:,:,:,grid%iold))

      do x = 1, nx
         do y = 1, ny

            rho_ = rho(y,x)
             ux_ =  ux(y,x)
             uy_ =  uy(y,x)

            pdf(y,x,:) = equilibrium(rho_, ux_, uy_)

         end do
      end do

      end associate

!      if (size(grid%f,4) > 2) then
!         grid%f(:,:,:,grid%imid) = grid%f(:,:,:,grid%iold)
!         grid%f(:,:,:,grid%inew) = grid%f(:,:,:,grid%iold)
!         call grid%collision()
!      end if

   end subroutine

   subroutine perform_step(grid)
      type(lattice_grid), intent(inout) :: grid

      call grid%streaming()  ! write from iold to inew
      call grid%collision()  ! update inew in place

      swap: block
         integer :: itmp
         itmp = grid%iold
         grid%iold = grid%inew
         grid%inew = itmp         
      end block swap

   end subroutine

   subroutine perform_triple_step(grid)
      type(lattice_grid), intent(inout) :: grid

      call grid%streaming()  ! write from iold to inew

      ! Store pre-collision PDF's
      grid%f(:,:,:,grid%iold) = grid%f(:,:,:,grid%inew)

      call grid%collision()  ! update inew in place

      swap: block
         integer :: itmp
         itmp = grid%iold
         grid%iold = grid%inew
         grid%inew = grid%imid
         grid%imid = itmp     
      end block swap

   end subroutine


   subroutine update_macros(grid)
      type(lattice_grid), intent(inout) :: grid

      integer :: ld

      ld = size(grid%f,1)

      call update_macros_kernel(grid%nx, grid%ny, ld, &
         grid%f(:,:,:,grid%inew), &
         grid%rho, grid%ux, grid%uy)

   contains

      subroutine update_macros_kernel(nx,ny,ld,f,grho,gux,guy)
         integer, intent(in) :: nx, ny, ld
         real(wp), intent(in) :: f(ld, nx, 0:8)
         real(wp), intent(inout), dimension(ny,nx) :: grho, gux, guy

         real(wp) :: rho, invrho, ux, uy, fs(0:8)
         integer :: x, y

         !$omp parallel do collapse(2) default(private) shared(nx,ny,f,grho,gux,guy)
         do x = 1, nx
            do y = 1, ny

               fs = f(y,x,:)

               ! density
               rho = fs(0) + (((fs(5) + fs(7)) + (fs(6) + fs(8))) + &
                      ((fs(1) + fs(3)) + (fs(2) + fs(4)))) 

               grho(y,x) = rho

               ! velocity
               invrho = 1.0_wp/rho
               ux = invrho * (((fs(5) - fs(7)) + (fs(8) - fs(6))) + (fs(1) - fs(3)))
               uy = invrho * (((fs(5) - fs(7)) + (fs(6) - fs(8))) + (fs(2) - fs(4)))

               gux(y,x) = ux
               guy(y,x) = uy

            end do
         end do
         !$omp end parallel do

      end subroutine

   end subroutine

   subroutine collide_bgk(grid)
      class(lattice_grid), intent(inout) :: grid
      integer :: ld

      ld = size(grid%f,1)

#if SPLIT
      call bgk_kernel_cache(grid%nx, grid%ny, ld, &
         grid%f(:,:,:,grid%inew), &
         grid%omega)
#else
      call bgk_kernel(grid%nx, grid%ny, ld, &
         grid%f(:,:,:,grid%inew), &
         grid%omega)
#endif

   contains

      subroutine bgk_kernel(nx,ny,ld,f1,omega)
         integer, intent(in) :: nx, ny, ld
         real(wp), intent(inout) :: f1(ld,nx,0:8)
         real(wp), intent(in) :: omega

         real(wp) :: rho, invrho, ux, uy
         real(wp) :: omegabar
         real(wp) :: fs(0:8), feq(0:8)

         integer :: x, y

         omegabar = 1.0_wp - omega

      !$omp parallel default(private) shared(f1,omega,omegabar,nx,ny)

         !$omp do collapse(2) schedule(static)
         do x = 1, nx
            do y = 1, ny

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
         !$omp end do

      !$omp end parallel

      end subroutine bgk_kernel

      subroutine bgk_kernel_cache(nx,ny,ld,f1,omega)
         integer, intent(in) :: nx,ny,ld
         real(wp), intent(inout) :: f1(ld,nx,0:8)
         real(wp), intent(in) :: omega

         real(wp) :: rho(ny), invrho, ux(ny), uy(ny), indp(ny)
         real(wp) :: fs(0:8)
         real(wp) :: omegabar, omega_w0, omega_ws, omega_wd

         real(wp), parameter :: one_third = 1.0_wp / 3.0_wp

         real(wp) :: vel_trm_13, vel_trm_24
         real(wp) :: vel_trm_57, vel_trm_68
         real(wp) :: velxpy, velxmy

         integer :: x, y


      !$omp parallel default(private) shared(f1,omega,nx,ny)
         
         omegabar = 1.0_wp - omega

         omega_w0 = 3.0_wp * omega * w0
         omega_ws = 3.0_wp * omega * ws
         omega_wd = 3.0_wp * omega * wd
         
         !$omp do schedule(static)
         do x = 1, nx

            do y = 1, ny
               ! pull pdfs travelling in different directions
               fs = f1(y,x,:)

               ! density
               rho(y) = (((fs(5) + fs(7)) + (fs(6) + fs(8))) + &
                         ((fs(1) + fs(3)) + (fs(2) + fs(4)))) + fs(0)

               invrho = 1.0_wp/rho(y)

               ! velocity
               ux(y) = invrho * (((fs(5) - fs(7)) + (fs(8) - fs(6))) + (fs(1) - fs(3)))
               uy(y) = invrho * (((fs(5) - fs(7)) + (fs(6) - fs(8))) + (fs(2) - fs(4)))

               indp(y) = one_third - 0.5_wp * (ux(y)**2 + uy(y)**2)

               ! update direction 0
               f1(y,x,0) = omegabar*fs(0) + omega_w0*rho(y)*indp(y)
            end do

            do y = 1, ny
            
               vel_trm_13 = indp(y) + 1.5_wp * ux(y) * ux(y)

               f1(y,x,1) = omegabar*f1(y,x,1) + omega_ws * rho(y) * (vel_trm_13 + ux(y))
               f1(y,x,3) = omegabar*f1(y,x,3) + omega_ws * rho(y) * (vel_trm_13 - ux(y))
            
            end do

            do y = 1, ny
               
               vel_trm_24 = indp(y) + 1.5_wp * uy(y) * uy(y)

               f1(y,x,2) = omegabar*f1(y,x,2) + omega_ws * rho(y) * (vel_trm_24 + uy(y))
               f1(y,x,4) = omegabar*f1(y,x,4) + omega_ws * rho(y) * (vel_trm_24 - uy(y))
            
            end do

            do y = 1, ny

               velxpy = ux(y) + uy(y)
               vel_trm_57 = indp(y) + 1.5_wp * velxpy * velxpy

               f1(y,x,5) = omegabar*f1(y,x,5) + omega_wd * rho(y) * (vel_trm_57 + velxpy)
               f1(y,x,7) = omegabar*f1(y,x,7) + omega_wd * rho(y) * (vel_trm_57 - velxpy)
            
            end do

            do y = 1, ny
               
               velxmy = ux(y) - uy(y)
               vel_trm_68 = indp(y) + 1.5_wp * velxmy * velxmy
               
               f1(y,x,6) = omegabar*f1(y,x,6) + omega_wd * rho(y) * (vel_trm_68 - velxmy)
               f1(y,x,8) = omegabar*f1(y,x,8) + omega_wd * rho(y) * (vel_trm_68 + velxmy)
            
            end do

         end do
         !$omp end do

      !$omp end parallel

      end subroutine bgk_kernel_cache

   end subroutine collide_bgk


   subroutine stream_fvm_bardow(grid)

      use interp, only: interp_hbox6, interp_vbox6

      class(lattice_grid), intent(inout) :: grid
      integer :: ld

      ld = size(grid%f,1)

      call fvm_bardow_kernel( &
         grid%nx, grid%ny, ld, &
         grid%f(:,:,:,grid%iold), &
         grid%f(:,:,:,grid%inew), &
         grid%dt)

   contains

      subroutine fvm_bardow_kernel(nx,ny,ld,fold,fnew,dt)
         integer, intent(in) :: nx, ny, ld
         real(wp), intent(in) :: fold(ld,nx,0:8)
         real(wp), intent(inout) :: fnew(ld,nx,0:8)
         real(wp), intent(in) :: dt

         integer :: q, x, y
         integer :: xp1, yp1, xm1, ym1

         real(wp) :: cxq, cyq
         real(wp) :: fc, fe, fn, fw, fs, fne, fnw, fsw, fse
         real(wp) :: cfe, cfn, cfw, cfs

         real(wp), parameter :: p2 = 0.5_wp, p8 = 0.125_wp, p6 = 1._wp/6._wp

      !$omp parallel default(private) shared(nx,ny,fold,fnew,dt)
         
         ! copy rest populations
         !$omp workshare
         fnew(1:ny,:,0) = fold(1:ny,:,0)
         !$omp end workshare

         do q = 1, 8
         
            cxq = dt*cx(q)
            cyq = dt*cy(q)

         !$omp do schedule(static)
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
         !$omp end do

         end do ! q
      
      !$omp end parallel

      end subroutine

   end subroutine

   subroutine stream_fdm_bardow(grid)
      class(lattice_grid), intent(inout) :: grid

      integer :: ld 

      ld = size(grid%f,1)

      call fdm_bardow_kernel( &
         grid%nx, grid%ny, ld, &
         grid%f(:,:,:,grid%iold), &
         grid%f(:,:,:,grid%inew), &
         grid%dt)

   contains

      subroutine fdm_bardow_kernel(nx,ny,ld,fold,fnew,dt)
         integer, intent(in) :: nx, ny,ld
         real(wp), intent(in) :: fold(ld,nx,0:8)
         real(wp), intent(inout) :: fnew(ld,nx,0:8)
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

         ! wls parameters
         real(wp), parameter :: one_third = 1._wp/3._wp


      !$omp parallel default(private) shared(nx,ny,fold,fnew,dt)

         ! copy rest populations
         !$omp workshare
         fnew(:,:,0) = fold(:,:,0)
         !$omp end workshare

         do q = 1, 8
         
            cxq = dt*cx(q)
            cyq = dt*cy(q)

            ! dt**2 is carried here
            cxxq = 0.5_wp*cxq*cxq
            cyyq = 0.5_wp*cyq*cyq
            cxyq = cxq*cyq

            !$omp do schedule(static)
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

#ifdef FDM_WLS
                  !> weighted least squares
                  !
                  !  no weight function
                  !  appears to be unstable
                  !  
                  dfx = one_sixth*((fne - fnw) + (fe - fw) + (fse - fsw))
                  dfy = one_sixth*((fne - fse) + (fn - fs) + (fnw - fsw))

                  dfxx = one_third*(fne - 2*fn + fnw) + one_third*(fe - 2*fc + fw) + one_third*(fse - 2*fs + fsw)
                  dfyy = one_third*(fne - 2*fe + fse) + one_third*(fn - 2*fc + fs) + one_third*(fnw - 2*fw + fsw)
                  dfxy = 0.25_wp*(fne - fnw + fsw - fse)
#elif FDM_WLS_GAUSS_V1
                  !> weighted least squares
                  !
                  !  gaussian weight, scaled to farthest node
                  !  
                  block
                     real(wp), parameter :: p1s = 0.2880584423829145035434_wp, &
                                            p1d = 0.1059707788085427065949_wp

                     real(wp), parameter :: p2c = -1.152233769531658458263_wp, &
                                            p2d1 = 0.5761168847658292291314_wp, &
                                            p2d2 = -0.4238831152341712149578_wp, &
                                            p2d = 0.2119415576170855242122_wp

                     dfx = p1s*(fe - fw) + p1d*(fne - fnw) + p1d*(fse - fsw)
                     dfy = p1s*(fn - fs) + p1d*(fne - fse) + p1d*(fnw - fsw)

                     dfxx = p2c*fc + p2d1*(fe + fw) + p2d2*(fn + fs) + p2d*(fne + fnw + fsw + fse)
                     dfyy = p2c*fc + p2d2*(fe + fw) + p2d1*(fn + fs) + p2d*(fne + fnw + fsw + fse)
                     dfxy = 0.25_wp*(fne - fnw + fsw - fse)
                  end block
#elif FDM_WLS_GAUSS_V2
                  !> weighted least squares
                  !
                  !  gaussian weight, scaled to closest node
                  !  match the values in the older rbff code
                  !  
                  block
                     real(wp), parameter :: p1s = 0.3934930210807994210853_wp, &
                                            p1d = 0.05325348945960039354075_wp

                     real(wp), parameter :: p2c = -1.573972084323197018207_wp, &
                                            p2d1 = 0.7869860421615988421706_wp, &
                                            p2d2 = -0.2130139578384016019186_wp, &
                                            p2d = 0.1065069789192007732037_wp

                     dfx = p1s*(fe - fw) + p1d*(fne - fnw) + p1d*(fse - fsw)
                     dfy = p1s*(fn - fs) + p1d*(fne - fse) + p1d*(fnw - fsw)

                     dfxx = p2c*fc + p2d1*(fe + fw) + p2d2*(fn + fs) + p2d*(fne + fnw + fsw + fse)
                     dfyy = p2c*fc + p2d2*(fe + fw) + p2d1*(fn + fs) + p2d*(fne + fnw + fsw + fse)
                     dfxy = 0.25_wp*(fne - fnw + fsw - fse)
                  end block
#elif FDM_ISO
                  !> isotropic derivatives, taken from
                  ! 
                  !   Kumar, A. (2004). Isotropic finite-differences. Journal of
                  !   Computational Physics, 201(1), 109-118, doi:10.1016/j.jcp.2004.05.005
                  !
                  !   the mixed derivative is already isotropic    
                  !
                  dfx = p2*(one_sixth*(fne - fnw) + two_thirds*(fe - fw) + one_sixth*(fse - fsw))
                  dfy = p2*(one_sixth*(fne - fse) + two_thirds*(fn - fs) + one_sixth*(fnw - fsw))

                  dfxx = one_twelth*(fne - 2*fn + fnw) + five_sixths*(fe - 2*fc + fw) + one_twelth*(fse - 2*fs + fsw)
                  dfyy = one_twelth*(fne - 2*fe + fse) + five_sixths*(fn - 2*fc + fs) + one_twelth*(fnw - 2*fw + fsw)
                  dfxy = 0.25_wp*(fne - fse - fnw + fsw)
#else 
                  ! first derivatives
                  dfx = p2*(fe - fw)
                  dfy = p2*(fn - fs)

                  ! second derivatives
                  dfxx = fe - 2*fc + fw
                  dfyy = fn - 2*fc + fs
                  dfxy = 0.25_wp*(fne - fse - fnw + fsw)
#endif
                  !
                  ! perform streaming
                  !
                  fnew(y,x,q) = fc - cxq*dfx - cyq*dfy + &
                     (cxxq*dfxx + cxyq*dfxy + cyyq*dfyy)

               end do
            end do
            !$omp end do

         end do

         !$omp end parallel
      end subroutine

   end subroutine


   subroutine stream_fdm_sofonea(grid)
      class(lattice_grid), intent(inout) :: grid

      integer :: ld

      ld = size(grid%f,1)

      call fdm_sofonea_kernel( &
         grid%nx, grid%ny, ld, &
         grid%f(:,:,:,grid%iold), &
         grid%f(:,:,:,grid%inew), &
         grid%dt)

   contains

      subroutine fdm_sofonea_kernel(nx,ny,ld,fold,fnew,dt)
         integer, intent(in) :: nx, ny, ld
         real(wp), intent(in) :: fold(ld,nx,0:8)
         real(wp), intent(inout) :: fnew(ld,nx,0:8)
         real(wp), intent(in) :: dt

         integer :: x, y
         integer :: xp1, yp1, xm1, ym1

         real(wp) :: du1, du2, fu, fc, fd

         real(wp), parameter :: p2 = 0.5_wp/sqrt(2._wp)

      !$omp parallel default(private) shared(nx,ny,fold,fnew,dt)

         ! copy rest populations
         !$omp workshare
         fnew(:,:,0) = fold(:,:,0)
         !$omp end workshare

         ! Direction 1
         !
         !$omp do schedule(static)
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
         !$omp end do

         ! Direction 2
         !
         !$omp do schedule(static)
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
         !$omp end do

         ! Direction 3
         !
         !$omp do schedule(static)
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
         !$omp end do

         ! Direction 4
         !
         !$omp do schedule(static)
         do x = 1, ny
            xp1 = mod(x, nx) + 1
            xm1 = mod(nx + x - 2, nx) + 1
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
         !$omp end do

         ! Direction 5
         !
         !$omp do schedule(static)
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
         !$omp end do

         ! Direction 6
         !
         !$omp do schedule(static)
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
         !$omp end do

         ! Direction 7
         !
         !$omp do schedule(static)
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
         !$omp end do

         ! Direction 8
         !
         !$omp do schedule(static)
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
         !$omp end do

      !$omp end parallel

      end subroutine

   end subroutine

   subroutine output_gnuplot(grid, step)
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
      
      call output_gnuplot_grid(fullname, &
         grid%nx,grid%ny, &
         grid%rho,grid%ux,grid%uy)

   end subroutine


   subroutine output_vtk(grid, step, binary)
      type(lattice_grid), intent(in) :: grid
      integer, intent(in), optional :: step
      logical, intent(in), optional :: binary

      logical :: binary_
      character(len=64) :: istr
      character(len=:), allocatable :: fullname

      binary_ = .false.
      if (present(binary)) binary_ = binary

      istr = ''
      if (present(step)) then
         write(istr,'(I0.9)') step
      end if

      fullname = ''
      if (allocated(grid%foldername)) then
         fullname = trim(grid%foldername) // '/'
      end if
      fullname = fullname // grid%filename // trim(istr) // '.vtk'

      if (binary_) then
         return
      else
         !call output_vtk_grid_ascii(fullname,grid%nx,grid%ny, &
         !   grid%rho,grid%ux,grid%uy)

         call output_vtk_structuredPoints(fullname,grid%nx,grid%ny, &
            grid%rho,grid%ux,grid%uy)
      end if

   end subroutine

   subroutine set_output_folder(grid,foldername)
      class(lattice_grid), intent(inout) :: grid
      character(len=*) :: foldername

      integer :: istat

      call execute_command_line("mkdir -p "//trim(foldername), &
         exitstat=istat, wait=.true.)

      if (istat /= 0) then
         write(*,'(A)') "[set_output_folder] error creating directory "//grid%foldername
         error stop
      end if

      grid%foldername = foldername

   end subroutine

end module



module periodic_dugks

   use precision, only: wp
   use fvm_bardow, only: lattice_grid, cx, cy

   implicit none
   private

   public :: perform_dugks_step

   public :: dugks_stream
   public :: dugks_collide


   ! Parameters related to the D2Q9 lattice

   real(wp), parameter :: one_third = 1.0_wp / 3.0_wp

   real(wp), parameter :: w0 = 4._wp / 9._wp, &
                          ws = 1._wp / 9._wp, &
                          wd = 1._wp / 36._wp

contains

   subroutine perform_dugks_step(grid)
      type(lattice_grid), intent(inout) :: grid

      call grid%collision()  ! update iold in place
      call grid%streaming()

      block
         integer :: itmp
         itmp = grid%iold
         grid%iold = grid%inew
         grid%inew = itmp
      end block

   end subroutine

   pure function tau_over_dt(grid)
      type(lattice_grid), intent(in) :: grid
      real(wp) :: tau_over_dt
      tau_over_dt = grid%tau/grid%dt
   end function

   subroutine dugks_collide(grid)
      class(lattice_grid), intent(inout) :: grid

      real(wp) :: omega, tau_d
      integer :: ld

      ld = size(grid%f,1)
      tau_d = tau_over_dt(grid)

      call copy_field(grid%nx, ld, &
         fsrc = grid%f(:,:,:,grid%iold), &
         fdst = grid%f(:,:,:,grid%inew))
      
      ! Collision for a full time-step
      omega = 1.0_wp/(tau_d + 0.5_wp)

      call kernel_bgk(grid%nx, grid%ny, ld, &
         grid%f(:,:,:,grid%inew), &
         grid%omega)

   ! TODO: we can save one collision, see original paper by Guo

#if DUGKS
      ! Collision for a half time-step
      omega = 0.75_wp*omega
#endif
      call kernel_bgk(grid%nx, grid%ny, ld, &
         grid%f(:,:,:,grid%iold), &
         omega)


   end subroutine


   subroutine kernel_bgk(nx,ny,ld,f,omega)
      integer, intent(in) :: nx,ny,ld
      real(wp), intent(inout) :: f(ld,nx,0:8)
      real(wp), intent(in) :: omega

      real(wp) :: rho(ny), invrho, ux(ny), uy(ny), indp(ny)
      real(wp) :: fs(0:8)
      real(wp) :: omegabar, omega_w0, omega_ws, omega_wd

      real(wp) :: vel_trm_13, vel_trm_24
      real(wp) :: vel_trm_57, vel_trm_68
      real(wp) :: velxpy, velxmy

      integer :: x, y

      !$omp parallel default(private) shared(f,omega,nx,ny)
         
      omegabar = 1.0_wp - omega

      omega_w0 = 3.0_wp * omega * w0
      omega_ws = 3.0_wp * omega * ws
      omega_wd = 3.0_wp * omega * wd
      
      !$omp do schedule(static)
      do x = 1, nx

         do y = 1, ny
            ! pull pdfs travelling in different directions
            fs = f(y,x,:)

            ! density
            rho(y) = (((fs(5) + fs(7)) + (fs(6) + fs(8))) + &
                      ((fs(1) + fs(3)) + (fs(2) + fs(4)))) + fs(0)

            invrho = 1.0_wp/rho(y)

            ! velocity
            ux(y) = invrho * (((fs(5) - fs(7)) + (fs(8) - fs(6))) + (fs(1) - fs(3)))
            uy(y) = invrho * (((fs(5) - fs(7)) + (fs(6) - fs(8))) + (fs(2) - fs(4)))

            indp(y) = one_third - 0.5_wp * (ux(y)**2 + uy(y)**2)

            ! update direction 0
            f(y,x,0) = omegabar*fs(0) + omega_w0*rho(y)*indp(y)
         end do

         do y = 1, ny
         
            vel_trm_13 = indp(y) + 1.5_wp * ux(y) * ux(y)

            f(y,x,1) = omegabar*f(y,x,1) + omega_ws * rho(y) * (vel_trm_13 + ux(y))
            f(y,x,3) = omegabar*f(y,x,3) + omega_ws * rho(y) * (vel_trm_13 - ux(y))
         
         end do

         do y = 1, ny
            
            vel_trm_24 = indp(y) + 1.5_wp * uy(y) * uy(y)

            f(y,x,2) = omegabar*f(y,x,2) + omega_ws * rho(y) * (vel_trm_24 + uy(y))
            f(y,x,4) = omegabar*f(y,x,4) + omega_ws * rho(y) * (vel_trm_24 - uy(y))
         
         end do

         do y = 1, ny

            velxpy = ux(y) + uy(y)
            vel_trm_57 = indp(y) + 1.5_wp * velxpy * velxpy

            f(y,x,5) = omegabar*f(y,x,5) + omega_wd * rho(y) * (vel_trm_57 + velxpy)
            f(y,x,7) = omegabar*f(y,x,7) + omega_wd * rho(y) * (vel_trm_57 - velxpy)
         
         end do

         do y = 1, ny
            
            velxmy = ux(y) - uy(y)
            vel_trm_68 = indp(y) + 1.5_wp * velxmy * velxmy
            
            f(y,x,6) = omegabar*f(y,x,6) + omega_wd * rho(y) * (vel_trm_68 - velxmy)
            f(y,x,8) = omegabar*f(y,x,8) + omega_wd * rho(y) * (vel_trm_68 + velxmy)
         
         end do

      end do
      !$omp end do

      !$omp end parallel

   end subroutine


   subroutine dugks_stream(grid)
      class(lattice_grid), intent(inout) :: grid

      integer :: ld
      real(wp) :: tau_d, omega

      ld = size(grid%f,1)

      tau_d = tau_over_dt(grid)
      omega = 1.0_wp/(4.0_wp*tau_d + 1.0_wp)

      call kernel_stream(grid%nx,grid%ny, ld, &
         grid%f(:,:,:,grid%iold), &
         grid%f(:,:,:,grid%inew), &
         grid%dt, omega)

   end subroutine

   subroutine kernel_stream(nx,ny,ld,ft,fp,dt,omega)
      integer, intent(in) :: nx, ny, ld
      real(wp), intent(in) :: ft(ld,nx,0:8)     ! \bar{f}^{+,n}
      real(wp), intent(inout) :: fp(ld,nx,0:8)  ! \tilde{f}^{+,n}

      real(wp), intent(in) :: dt, omega

      real(wp) :: cfw(ny,0:8), cfn(ny,0:8), cfe(ny,0:8), cfs(ny,0:8)
      
#if DUGKS
      real(wp) :: cfrho(ny), cfux(ny), cfuy(ny), cfindp(ny)
#endif

      integer :: x, y, q
      integer :: xp1, xm1, yp1, ym1

      real(wp) :: cxq, cyq
      
      real(wp), parameter :: p2 = 0.5_wp, &
                             p8 = 0.125_wp
      
      real(wp) :: fc, fe, fn, fw, fs, fne, fnw, fsw, fse

      real(wp) :: fluxW, fluxN, fluxE, fluxS

      !$omp parallel default(private) shared(nx,ny,ld,ft,fp,dt,omega)
      !$omp do schedule(static)
      do x = 1, nx

         ! locate neighbor nodes
         xp1 = mod(x, nx) + 1
         xm1 = mod(nx + x - 2, nx) + 1

         ! Interpolate values at the WEST and NORTH cell faces
         ! 
         ! All directions are needed, so we can calculate the 
         ! macroscopic values
         !  
         do q = 0, 8
         
         cxq = dt*cx(q)
         cyq = dt*cy(q)

         do y = 1, ny

            yp1 = mod(y, ny) + 1
            ym1 = mod(ny + y - 2, ny) + 1

            ! read pdf values
            fc  = ft(y  , x  , q)
            fe  = ft(y  , xp1, q)
            fn  = ft(yp1, x  , q)
            fw  = ft(y  , xm1, q)
            fs  = ft(ym1, x  , q)
            fne = ft(yp1, xp1, q)
            fnw = ft(yp1, xm1, q)
            fsw = ft(ym1, xm1, q)
            fse = ft(ym1, xp1, q)

            ! WEST cell face
            cfw(y,q) = p2*(fc + fw) - p2*cxq*(fc - fw) - &
                       p8*cyq*(fnw + fn - fsw - fs)

            ! NORTH cell face
            cfn(y,q) = p2*(fc + fn) - p2*cyq*(fn - fc) - &
                       p8*cxq*(fne + fe - fnw - fw)

            ! EAST cell face
            cfe(y,q) = p2*(fc + fe) - p2*cxq*(fe - fc) - &
                       p8*cyq*(fne + fn - fse - fs)

            ! SOUTH cell face
            cfs(y,q) = p2*(fc + fs) - p2*cyq*(fc - fs) - &
                       p8*cxq*(fse + fe - fsw - fw)
         end do
         end do

#if DUGKS
         ! WEST face
         call update_ew(ny,cfw,cfrho,cfux,cfuy,cfindp,omega)
         ! EAST face 
         call update_ew(ny,cfe,cfrho,cfux,cfuy,cfindp,omega)
         ! NORTH face
         call update_ns(ny,cfn,cfrho,cfux,cfuy,cfindp,omega)
         ! SOUTH face
         call update_ns(ny,cfs,cfrho,cfux,cfuy,cfindp,omega)
#endif

         ! TODO: split this into two updates
         !       one for the east & west faces, and
         !       one for the north & south faces
         !
         do q = 1, 8
         
         cxq = dt*cx(q)
         cyq = dt*cy(q)

         do y = 1, ny

            yp1 = mod(y, ny) + 1
            ym1 = mod(ny + y - 2, ny) + 1

            fluxE = cfe(y,q)
            fluxW = cfw(y,q)
            fluxN = cfn(y,q)
            fluxS = cfs(y,q)

            fp(y,x,q) = fp(y,x,q) - cxq*(fluxE - fluxW) - cyq*(fluxN - fluxS)

         end do
         end do

      end do
      !$omp end do
      !$omp end parallel

   contains

#if DUGKS

      subroutine update_ew(ny,f,rho,ux,uy,indp,omega)
         integer, intent(in) :: ny
         real(wp), intent(inout) :: f(ny,0:8)
         real(wp), intent(out) :: rho(ny), ux(ny), uy(ny), indp(ny)
         real(wp), intent(in) :: omega

         integer :: y
         real(wp) :: invrho
         real(wp) :: omegabar, omega_w0, omega_ws, omega_wd

         real(wp) :: vel_trm_13, vel_trm_57, vel_trm_68
         real(wp) :: velxpy, velxmy

         omegabar = 1.0_wp - omega

         omega_w0 = 3.0_wp * omega * w0
         omega_ws = 3.0_wp * omega * ws
         omega_wd = 3.0_wp * omega * wd

         do y = 1, ny

            rho(y) = (((f(y,5) + f(y,7)) + (f(y,6) + f(y,8))) + &
                      ((f(y,1) + f(y,3)) + (f(y,2) + f(y,4)))) + f(y,0)
            invrho = 1.0_wp/rho(y)

            ux(y) = invrho * (((f(y,5) - f(y,7)) + (f(y,8) - f(y,6))) + (f(y,1) - f(y,3)))
            uy(y) = invrho * (((f(y,5) - f(y,7)) + (f(y,6) - f(y,8))) + (f(y,2) - f(y,4)))

            indp(y) = one_third - 0.5_wp * (ux(y)**2 + uy(y)**2)
         end do

         do y = 1, ny
         
            vel_trm_13 = indp(y) + 1.5_wp * ux(y) * ux(y)

            f(y,1) = omegabar*f(y,1) + omega_ws * rho(y) * (vel_trm_13 + ux(y))
            f(y,3) = omegabar*f(y,3) + omega_ws * rho(y) * (vel_trm_13 - ux(y))
         
         end do

         do y = 1, ny

            velxpy = ux(y) + uy(y)
            vel_trm_57 = indp(y) + 1.5_wp * velxpy * velxpy

            f(y,5) = omegabar*f(y,5) + omega_wd * rho(y) * (vel_trm_57 + velxpy)
            f(y,7) = omegabar*f(y,7) + omega_wd * rho(y) * (vel_trm_57 - velxpy)
         
         end do

         do y = 1, ny
            
            velxmy = ux(y) - uy(y)
            vel_trm_68 = indp(y) + 1.5_wp * velxmy * velxmy
            
            f(y,6) = omegabar*f(y,6) + omega_wd * rho(y) * (vel_trm_68 - velxmy)
            f(y,8) = omegabar*f(y,8) + omega_wd * rho(y) * (vel_trm_68 + velxmy)
         
         end do

      end subroutine

      subroutine update_ns(ny,f,rho,ux,uy,indp,omega)
         integer, intent(in) :: ny
         real(wp), intent(inout) :: f(ny,0:8)
         real(wp), intent(out) :: rho(ny), ux(ny), uy(ny), indp(ny)
         real(wp), intent(in) :: omega

         integer :: y
         real(wp) :: invrho
         real(wp) :: omegabar, omega_w0, omega_ws, omega_wd

         real(wp) :: vel_trm_24, vel_trm_57, vel_trm_68
         real(wp) :: velxpy, velxmy

         omegabar = 1.0_wp - omega

         omega_w0 = 3.0_wp * omega * w0
         omega_ws = 3.0_wp * omega * ws
         omega_wd = 3.0_wp * omega * wd

         do y = 1, ny

            rho(y) = (((f(y,5) + f(y,7)) + (f(y,6) + f(y,8))) + &
                      ((f(y,1) + f(y,3)) + (f(y,2) + f(y,4)))) + f(y,0)
            
            invrho = 1.0_wp/rho(y)

            ux(y) = invrho * (((f(y,5) - f(y,7)) + (f(y,8) - f(y,6))) + (f(y,1) - f(y,3)))
            uy(y) = invrho * (((f(y,5) - f(y,7)) + (f(y,6) - f(y,8))) + (f(y,2) - f(y,4)))

            indp(y) = one_third - 0.5_wp * (ux(y)**2 + uy(y)**2)
         
         end do

         do y = 1, ny
            
            vel_trm_24 = indp(y) + 1.5_wp * uy(y) * uy(y)

            f(y,2) = omegabar*f(y,2) + omega_ws * rho(y) * (vel_trm_24 + uy(y))
            f(y,4) = omegabar*f(y,4) + omega_ws * rho(y) * (vel_trm_24 - uy(y))
         
         end do

         do y = 1, ny

            velxpy = ux(y) + uy(y)
            vel_trm_57 = indp(y) + 1.5_wp * velxpy * velxpy

            f(y,5) = omegabar*f(y,5) + omega_wd * rho(y) * (vel_trm_57 + velxpy)
            f(y,7) = omegabar*f(y,7) + omega_wd * rho(y) * (vel_trm_57 - velxpy)
         
         end do

         do y = 1, ny
            
            velxmy = ux(y) - uy(y)
            vel_trm_68 = indp(y) + 1.5_wp * velxmy * velxmy
            
            f(y,6) = omegabar*f(y,6) + omega_wd * rho(y) * (vel_trm_68 - velxmy)
            f(y,8) = omegabar*f(y,8) + omega_wd * rho(y) * (vel_trm_68 + velxmy)
         
         end do

      end subroutine
      
#endif

   end subroutine


   subroutine copy_field(nx,ld,fsrc,fdst)
      integer, intent(in) :: nx, ld
      real(wp), intent(in) :: fsrc(ld,nx,0:8)
      real(wp), intent(out) :: fdst(ld,nx,0:8)

      !$omp parallel default(private) shared(fsrc,fdst)

      !$omp workshare
      fdst(:,:,0) = fsrc(:,:,0)
      !$omp end workshare

      !$omp workshare
      fdst(:,:,1) = fsrc(:,:,1)
      !$omp end workshare

      !$omp workshare
      fdst(:,:,2) = fsrc(:,:,2)
      !$omp end workshare

      !$omp workshare
      fdst(:,:,3) = fsrc(:,:,3)
      !$omp end workshare

      !$omp workshare
      fdst(:,:,4) = fsrc(:,:,4)
      !$omp end workshare

      !$omp workshare
      fdst(:,:,5) = fsrc(:,:,5)
      !$omp end workshare

      !$omp workshare
      fdst(:,:,6) = fsrc(:,:,6)
      !$omp end workshare

      !$omp workshare
      fdst(:,:,7) = fsrc(:,:,7)
      !$omp end workshare

      !$omp workshare
      fdst(:,:,8) = fsrc(:,:,8)
      !$omp end workshare

      !$omp end parallel

   end subroutine copy_field

end module
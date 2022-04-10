module collision_bgk

   use precision, only: wp
   use fvm_bardow, only: lattice_grid, equilibrium

   implicit none
   private

   public :: collide_bgk

   real(wp), parameter :: w0 = 4._wp / 9._wp, &
                          ws = 1._wp / 9._wp, &
                          wd = 1._wp / 36._wp
                          
contains

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

end module
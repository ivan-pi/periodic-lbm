module collision_bgk_improved

   use precision, only: wp
   use fvm_bardow, only: lattice_grid

   implicit none
   private

   public :: collide_bgk_improved

   real(wp), parameter :: one_third = 1.0_wp/3.0_wp
   real(wp), parameter :: two_thirds = 2.0_wp/3.0_wp

contains

   subroutine collide_bgk_improved(grid)
      class(lattice_grid), intent(inout) :: grid

      call bgk_improved_kernel(grid%nx, grid%ny,&
         grid%f(:,:,:,grid%inew), grid%omega)

   contains

      subroutine bgk_improved_kernel(nx,ny,f1,omega)
         integer, intent(in) :: nx, ny
         real(wp), intent(inout) :: f1(ny,nx,0:8)
         real(wp), intent(in) :: omega

         real(wp) :: fac, omegabar
         real(wp) :: vC, vE, vN, vW, vS, vNE, vNW, vSW, vSE

         real(wp) :: sumx1, sumxn, sumy1, sumyn
         real(wp) :: m10, m01, m20, m02, u2, v2, Gx, Gy
         real(wp) :: X0, X1, XN, Y0, Y1, YN
         real(wp) :: rho, invrho, rho_omega

         integer :: x, y

         fac = 4.5_wp - 2.25_wp * omega
         omegabar = 1.0_wp - omega

         do y = 1, ny
            do x = 1, nx

               vC  = f1(y,x,0)
               VE  = f1(y,x,1)
               vN  = f1(y,x,2)
               vW  = f1(y,x,3)
               vS  = f1(y,x,4)
               vNE = f1(y,x,5)
               vNW = f1(y,x,6)
               vSW = f1(y,x,7)
               vSE = f1(y,x,8)

               rho = (((vNE + vSW) + (vNW + vSE)) + &
                      ((vE + vW) + (vN + vS))) + vC

               invrho = 1.0_wp/rho

               !ux = invrho * (((vNE - vSW) + (vSE - vNW)) + (vE - vW))
               !uy = invrho * (((vNE - vSW) + (vNW - vSE)) + (vN - vS))

               sumX1 = vE + vNE + vSE
               sumXN = vW + vNW + vSW

               sumY1 = VN + vNE + vNW
               sumYN = vS + vSE + vSW

               m10 = invrho * (sumX1 - sumXN)
               m01 = invrho * (sumY1 - sumYN)

               u2 = m10*m10
               v2 = m01*m01

               m20 = invrho * (sumX1 + sumXN)
               m02 = invrho * (sumY1 + sumYN)

               Gx = fac * u2 * (m20 - one_third - u2)
               Gy = fac * v2 * (m02 - one_third - v2)

               X0 = -two_thirds + u2 + Gx
               X1 = - (X0 + 1.0_wp + m10) * 0.5_wp
               XN = X1 + m10

               Y0 = -two_thirds + v2 + Gy
               Y1 = - (Y0 + 1.0_wp + m01) * 0.5_wp
               YN = Y1 + m01

               rho_omega = rho * omega
               X0 = X0 * rho_omega
               X1 = X1 * rho_omega
               XN = XN * rho_omega

               f1(y,x,0) = omegabar * vC  + X0*Y0
               f1(y,x,1) = omegabar * vE  + X1*Y0
               f1(y,x,2) = omegabar * vN  + X0*Y1
               f1(y,x,3) = omegabar * vW  + XN*Y0
               f1(y,x,4) = omegabar * vS  + X0*YN
               f1(y,x,5) = omegabar * vNE + X1*Y1
               f1(y,x,6) = omegabar * vNW + XN*Y1
               f1(y,x,7) = omegabar * vSW + XN*YN
               f1(y,x,8) = omegabar * vSE + X1*YN

            end do
         end do

      end subroutine

   end subroutine

end module
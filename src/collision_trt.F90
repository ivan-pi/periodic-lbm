module collision_trt

   use precision, only: wp
   use fvm_bardow, only: lattice_grid

   implicit none
   private

   public :: magic_number, lambda_d
   public :: collide_trt

   ! common prefactors for calculating the equilibrium parts
   real(wp), parameter :: t0 = 4.0_wp / 9.0_wp
   real(wp), parameter :: t1x2 = ( 1.0_wp / 9.0_wp  ) * 2.0_wp
   real(wp), parameter :: t2x2 = ( 1.0_wp / 36.0_wp ) * 2.0_wp

   ! speed of sound related factors
   real(wp), parameter :: inv2csq2 = 1.0_wp / ( 2.0_wp * ( 1.0_wp / 3.0_wp ) * ( 1.0_wp / 3.0_wp ) )
   real(wp), parameter :: fac1 = t1x2 * inv2csq2
   real(wp), parameter :: fac2 = t2x2 * inv2csq2

contains

   pure real(wp) function magic_number(le,ld)
      !> Dimensionless relaxation rates
      real(wp), intent(in) :: le, ld
      magic_number = (2.0_wp - le) * (2.0_wp - ld) / (4.0_wp*le*ld)
   end function

   pure real(wp) function lambda_d(omega,x)
      real(wp), intent(in) :: omega, x
      lambda_d = (4.0_wp - 2.0_wp*omega) / &
         (4.0_wp * x * omega + 2.0_wp - omega)
   end function


   ! TRT Kernel
   !
   ! Incompressible equilibrium
   !
   subroutine collide_trt(grid)
      class(lattice_grid), intent(inout) :: grid

      integer :: ld
      real(wp) :: lambda_even, lambda_odd

      lambda_even = grid%omega
      lambda_odd  = lambda_d(grid%omega, grid%trt_magic)

      ld = size(grid%f,1)

#if SPLIT
      call trt_split(grid%nx, grid%ny, &
         grid%f(:,:,:,grid%inew), ld, &
         lambda_even, lambda_odd)
#else
      call trt_naive(grid%nx, grid%ny, &
         grid%f(:,:,:,grid%inew), ld, &
         lambda_even, lambda_odd)
#endif

   contains

      subroutine trt_naive(nx,ny,f1,ldf1,lambda_e,lambda_d)
         integer, intent(in) :: nx, ny, ldf1
         real(wp), intent(inout) :: f1(ldf1,nx,0:8)
         real(wp), intent(in) :: lambda_e, lambda_d

         
         real(wp) :: rho, invrho, velX, velY

         ! relaxation parameter variables
         real(wp) :: lambda_e_scaled, lambda_d_scaled
         
         real(wp) :: sym_NE_SW, asym_NE_SW
         real(wp) :: sym_SE_NW, asym_SE_NW

         real(wp) :: sym_N_S, asym_N_S
         real(wp) :: sym_E_W, asym_E_W

         real(wp) :: velX2, velY2
         real(wp) :: velXPY, velXMY
         real(wp) :: vC, vE, vN, vW, vS, vNE, vNW, vSW, vSE
         real(wp) :: feq_common

         integer :: x, y

         !$omp parallel default(private) shared(f1,lambda_e,lambda_d,nx,ny)

         lambda_e_scaled = 0.5_wp * lambda_e   ! 0.5 times the usual value
         lambda_d_scaled = 0.5_wp * lambda_d   ! ... due to the way of calculations

         !$omp do schedule(static)
         do x = 1, nx  
            do y = 1, ny

            ! pull pdfs
            vC  = f1(y,x,0)
            VE  = f1(y,x,1)
            vN  = f1(y,x,2)
            vW  = f1(y,x,3)
            vS  = f1(y,x,4)
            vNE = f1(y,x,5)
            vNW = f1(y,x,6)
            vSW = f1(y,x,7)
            vSE = f1(y,x,8)

            !
            ! macroscopic values (incompressible equilibrium)
            !

            rho = (((vNE + vSW) + (vNW + vSE)) + &
                   ((vE + vW) + (vN + vS))) + vC

            velX = (((vNE - vSW) + (vSE - vNW)) + (vE - vW))
            velY = (((vNE - vSW) + (vNW - vSE)) + (vN - vS))

            velX2 = velX*velX
            velY2 = velY*velY
            
            !
            ! compute new populations
            !

            feq_common = rho - 1.5_wp*(velX2 + velY2)

            f1(y,x,0) = vC * (1.0_wp - lambda_e) + lambda_e * t0 * feq_common

            velXPY = velX + velY

             sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common )
            asym_NE_SW = lambda_d_scaled * ( vNE - vSW - 3.0_wp * t2x2 * velXPY )
            f1(y,x,5) = vNE - sym_NE_SW - asym_NE_SW
            f1(y,x,7) = vSW - sym_NE_SW + asym_NE_SW

            velXMY = velX - velY

             sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common )
            asym_SE_NW = lambda_d_scaled * ( vSE - vNW - 3.0_wp * t2x2 * velXMY ) 
            f1(y,x,8) = vSE - sym_SE_NW - asym_SE_NW
            f1(y,x,6) = vNW - sym_SE_NW + asym_SE_NW

             sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY2 - t1x2 * feq_common )
            asym_N_S = lambda_d_scaled * ( vN - vS - 3.0_wp * t1x2 * velY )
            f1(y,x,2) = vN - sym_N_S - asym_N_S
            f1(y,x,4) = vS - sym_N_S + asym_N_S

             sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX2 - t1x2 * feq_common )
            asym_E_W = lambda_d_scaled * ( vE - vW - 3.0_wp * t1x2 * velX )

            f1(y,x,1) = vE - sym_E_W - asym_E_W
            f1(y,x,3) = vW - sym_E_W + asym_E_W

            end do
         end do
         !$omp end do

         !$omp end parallel

      end subroutine

      subroutine trt_split(nx,ny,f1,ldf1,lambda_e,lambda_d)
         integer, intent(in) :: nx, ny, ldf1
         real(wp), intent(inout) :: f1(ldf1,nx,0:8)

         real(wp), intent(in) :: lambda_e, lambda_d

         real(wp) :: rho(ny), velX(ny), velY(ny)
         !dir$ attributes align:64 :: rho
         !dir$ attributes align:64 :: velX
         !dir$ attributes align:64 :: velY

         ! relaxation parameter variables
         real(wp) :: lambda_e_scaled, lambda_d_scaled
         
         real(wp) :: sym_NE_SW, asym_NE_SW
         real(wp) :: sym_SE_NW, asym_SE_NW

         real(wp) :: sym_N_S, asym_N_S
         real(wp) :: sym_E_W, asym_E_W

         real(wp) :: velXPY, velXMY
         real(wp) :: vC, vE, vN, vW, vS, vNE, vNW, vSW, vSE
         
         real(wp) :: feq_common(ny)
         !dir$ attributes align: 64 :: feq_common

         integer :: x, y

         !dir$ assume_aligned f1: 64

      !$omp parallel default(private) shared(f1,lambda_e,lambda_d,nx,ny)

         lambda_e_scaled = 0.5_wp * lambda_e   ! 0.5 times the usual value
         lambda_d_scaled = 0.5_wp * lambda_d   ! ... due to the way of calculations

         !$omp do schedule(static)
         do x = 1, nx
            do y = 1, ny

               ! pull pdfs
               vC  = f1(y,x,0)
               VE  = f1(y,x,1)
               vN  = f1(y,x,2)
               vW  = f1(y,x,3)
               vS  = f1(y,x,4)
               vNE = f1(y,x,5)
               vNW = f1(y,x,6)
               vSW = f1(y,x,7)
               vSE = f1(y,x,8)

               !
               ! macroscopic values (incompressible equilibrium)
               !

               rho(y) = (((vNE + vSW) + (vNW + vSE)) + &
                      ((vE + vW) + (vN + vS))) + vC

               velX(y) = (((vNE - vSW) + (vSE - vNW)) + (vE - vW))
               velY(y) = (((vNE - vSW) + (vNW - vSE)) + (vN - vS))
               !
               ! compute new populations
               !

               feq_common(y) = rho(y) - 1.5_wp*(velX(y)*velX(y) + velY(y)*velY(y))
               f1(y,x,0) = vC * (1.0_wp - lambda_e) + lambda_e * t0 * feq_common(y)
            
            end do

            do y = 1, ny

               vNE = f1(y,x,5)
               vSW = f1(y,x,7)
         
               velXPY = velX(y) + velY(y)

               sym_NE_SW  = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common(y) )
               asym_NE_SW = lambda_d_scaled * ( vNE - vSW - 3.0_wp * t2x2 * velXPY )

               f1(y,x,5) = f1(y,x,5) - sym_NE_SW - asym_NE_SW
               f1(y,x,7) = f1(y,x,7) - sym_NE_SW + asym_NE_SW

            end do

            do y = 1, ny
               
               vNW = f1(y,x,6)
               vSE = f1(y,x,8)

               velXMY = velX(y) - velY(y)
               
                sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common(y) )
               asym_SE_NW = lambda_d_scaled * ( vSE - vNW - 3.0_wp * t2x2 * velXMY ) 
               
               f1(y,x,8) = f1(y,x,8) - sym_SE_NW - asym_SE_NW
               f1(y,x,6) = f1(y,x,6) - sym_SE_NW + asym_SE_NW

            end do

            do y = 1, ny

               vN = f1(y,x,2)
               vS = f1(y,x,4)
               
                sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY(y) * velY(y) - t1x2 * feq_common(y) )
               asym_N_S = lambda_d_scaled * ( vN - vS - 3.0_wp * t1x2 * velY(y) )
               
               f1(y,x,2) = f1(y,x,2) - sym_N_S - asym_N_S
               f1(y,x,4) = f1(y,x,4) - sym_N_S + asym_N_S

            end do

            do y = 1, ny

               vE = f1(y,x,1)
               vW = f1(y,x,3)
               
                sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX(y) * velX(y) - t1x2 * feq_common(y) )
               asym_E_W = lambda_d_scaled * ( vE - vW - 3.0_wp * t1x2 * velX(y) )
               
               f1(y,x,1) = f1(y,x,1) - sym_E_W - asym_E_W
               f1(y,x,3) = f1(y,x,3) - sym_E_W + asym_E_W

            end do
         end do
         !$omp end do

      !$omp end parallel

      end subroutine


   end subroutine

end module

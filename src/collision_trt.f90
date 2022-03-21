module collision_trt

   implicit none


contains

   subroutine collide_trt(grid)
      class(lattice_grid), intent(inout) :: grid

      real(wp) :: lambda_e, lambda_d



      call trt_kernel(grid%nx, grid%ny,&
         grid%f(:,:,:,grid%inew), &
         lambda_e, lambda_d)

   contains

      subroutine trt_kernel(nx,ny,f1,omega)
         integer, intent(in) :: nx, ny
         real(wp), intent(inout) :: f1(ny,nx,0:8)
         real(wp), intent(in) :: omega

         real(wp) :: rho,velX, velY
         real(wp) :: omegabar
         real(wp) :: fs(0:8), feq(0:8)

         integer :: x, y

         omegabar = 1.0_wp - omega

         !$omp parallel do collapse(2) default(private) shared(f1,omega,omegabar)
         do y = 1, ny
            do x = 1, nx  

  ! common prefactors for calculating the equilibrium parts
  real(wp), parameter :: t0 = 4.0_wp/9.0_wp
  real(wp), parameter :: t1x2 = 1.0_wp/9.0_wp * 2.0_wp
  real(wp), parameter :: t2x2 = 1.0_wp/36.0_wp * 2.0_wp

  ! speed of sound related factors
  real(wp), parameter :: inv2csq2 = 1.0_wp / ( 2.0_wp * ( 1.0_wp / 3.0_wp ) * ( 1.0_wp / 3.0_wp ) )
  real(wp), parameter :: fac1 = t1x2 * inv2csq2
  real(wp), parameter :: fac2 = t2x2 * inv2csq2

  ! relaxation parameter variables
  real(wp) :: lambda_e_scaled, lambda_d_scaled
  
  real(wp) :: sym_NE_SW, asym_NE_SW
  real(wp) :: sym_SE_NW, asym_SE_NW
  
  real(wp) :: sym_N_S, asym_N_S
  real(wp) :: sym_E_W, asym_E_W

  real(wp) :: velXPY, velXMY

  lambda_e_scaled = 0.5_wp * lambda_e   ! 0.5 times the usual value
  lambda_d_scaled = 0.5_wp * lambda_d   ! ... due to the way of calculations

  feq_common = rho - 1.5_wp*(ux*ux + uy*uy)

  dst(0) = vC * (1.0_wp - lambda_e) + lambda_e * t0 * feq_common

  velXPY = velX + velY
   
   sym_NE_SW = lambda_e_scaled * ( vNE + vSW - fac2 * velXPY * velXPY - t2x2 * feq_common )
  asym_NE_SW = lambda_d_scaled * ( vNE - vSW - 3.0_wp * t2x2 * velXPY )
  dst(5) = vNE - sym_NE_SW - asym_NE_SW
  dst(7) = vSW - sym_NE_SW + asym_NE_SW

  velXMY = velX - velY

   sym_SE_NW = lambda_e_scaled * ( vSE + vNW - fac2 * velXMY * velXMY - t2x2 * feq_common )
  asym_SE_NW = lambda_d_scaled * ( vSE - vNW - 3.0_wp * t2x2 * velXMY ) 
  dst(8) = vSE - sym_SE_NW - asym_SE_NW
  dst(6) = vNW - sym_SE_NW + asym_SE_MW

   sym_N_S = lambda_e_scaled * ( vN + vS - fac1 * velY * velY - t1x2 * feq_common )
  asym_N_S = lambda_d_scaled * ( vN - vS - 3.0_wp * t1x2 * velY )
  dst(2) = vN - sym_N_S - asym_N_S
  dst(4) = vS - sym_N_S + asym_N_S

   sym_E_W = lambda_e_scaled * ( vE + vW - fac1 * velX * velX - t1x2 * feq_common )
  asym_E_W = lambda_d_scaled * ( vE - vW - 3.0_wp * t1x2 * velX )

  dst(1) = vE - sym_E_W - asym_E_W
  dst(3) = vW - sym_E_W + asym_E_W

end module

real(wp), parameter :: one_third = 1.0_wp / 3.0_wp

real(wp) :: omega_trm, omega_w0, omega_w1, omega_w2

real(wp) :: velXX, velYY, velXmY, velXpY

real(wp) :: dir_indep_trm, &
            vel_trm_E_W, vel_trm_N_S, &
            vel_trm_NW_SE, vel_trm_NE_SW

omega_trm = 1.0_wp - omega
omega_w0  = 3.0_wp * w0 * omega
omega_w1  = 3.0_wp * w1 * omega
omega_w2  = 3.0_wp * w2 * omega

! bgk

velXX = velX * velX
velYY = velY * velY

dir_indep_trm = one_third * rho - 0.5_wp * (velXX + velYY)

dst(0) = omega_trm * vC + omega_w0 * dir_indep_trm

vel_trm_E_W = dir_indep_trm + 1.5_wp * velXX
vel_trm_N_S = dir_indep_trm + 1.5_wp * velYY

dst(1) = omega_trm * vE + omega_w1 * ( vel_trm_E_W + velX)
dst(3) = omega_trm * vW + omega_w1 * ( vel_trm_E_W - velX)
dst(2) = omega_trm * vN + omega_w1 * ( vel_trm_N_S + velY)
dst(4) = omega_trm * vS + omega_w1 * ( vel_trm_N_S - velY)

velXmY = velX - velY
vel_trm_NW_SE = dir_indep_trm + 1.5_wp * velXmY * velXmY

dst(6) = omega_trm * vNW + omega_w2 * ( vel_trm_NW_SE - velXmY )
dst(8) = omega_trm * vSE + omega_w2 * ( vel_trm_NW_SE + velXmY )

velXpY = velX + velY
vel_trm_NE_SW = dir_indep_trm + 1.5_wp * velXpY * velXpY

dst(5) = omega_trm * vNE + omega_w2 * ( vel_trm_NE_SW + velXpY )
dst(7) = omega_trm * vNE + omega_w2 * ( vel_trm_NE_SW - velXpY )
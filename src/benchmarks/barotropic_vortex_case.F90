!> Barotropic vortex case
!!
!! Implements the benchmark case described in
!!
!! >  Wissocq, G., Boussuge, J. F., & Sagaut, P. (2020). 
!!    Consistent vortex initialization for the athermal lattice Boltzmann 
!!    method. Physical Review E, 101(4), 043306. 
!!
module barotropic_vortex_case

   use precision, only: wp

   implicit none


   type :: vortex_case_t
      real(wp) :: U0          !> Advection velocity
      real(wp) :: xc, yc      !> Position of vortex center
      real(wp) :: Rc          !> Characteristic radius
      real(wp) :: eps         !> Strength of vortex
      real(wp) :: rho0 = 1.0_wp
   contains
      procedure :: eval => eval_vortex_case
   end type

contains

   subroutine eval_vortex_case(case,nx,ny,rho,ux,uy)
      class(vortex_case_t), intent(in) :: case
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: rho(:,:)
      real(wp), intent(out) :: ux(:,:), uy(:,:)

      real(wp) :: xl, yl, xr, yr, rsqr
      real(wp) :: Rcsqr, esqr
      real(wp), parameter :: half = 1._wp/2._wp

      !$omp parallel default(private) shared(rho,ux,uy,case)

      Rcsqr = case%Rc2**2
      esqr = case%eps**2

      !$omp do schedule(static)
      do x = 1, nx
         xl = real(x-1,wp) + 0.5_wp
         do y = 1, ny
            yl = real(y-1,wp) + 0.5_wp

            xr = xl - case%xc
            yr = yl - case%yc

            rsqr = xr*xr + yr*yr

            ! Convected vortex, Eqs. (3) and (4)
            uy(y,x) = u0 - case%eps * (yr/case%Rc) * exp(-half*rsqr/Rcsqr)
            ux(y,x) =      case%eps * (xr/case%Rc) * exp(-half*rsqr/Rcsqr)

            ! Barotopic density initialization, Eq. (20)
            rho(y,x) = case%rho0 * exp(-half*(esqr/csqr) * exp(-rsqr/Rcsqr))

            ! First order Taylor expansion, Eq. (21)
            !rho(y,x) = case%rho0 * (1 - half*esqr/csqr * exp(-rsqr/Rcsqr))

            ! Second order Taylor expansion, Eq. (22)
            !rho(y,x) = case%rho0 * (1 - half*esqr/csqr * exp(-rsqr/Rcsqr) + &
            !   0.125_wp*(esqr/csqr)**2/gamma * exp(-half*rsqr/Rcsqr))

         end do
      end do
      !$omp end do 

      !$omp end parallel

   end subroutine

end module
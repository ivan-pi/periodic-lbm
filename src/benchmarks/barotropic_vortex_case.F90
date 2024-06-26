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


   ! TODO: Think about whether csqr should be part
   !       of the derived type, or the eval procedure

   type :: vortex_case_t
      real(wp) :: U0          !> Advection velocity
      real(wp) :: xc, yc      !> Position of vortex center
      real(wp) :: Rc          !> Characteristic radius
      real(wp) :: eps         !> Strength of vortex
      real(wp) :: rho0 = 1.0_wp
      real(wp) :: csqr = 1._wp/3._wp
   contains
      procedure :: eval => eval_vortex_case, &
                           eval_vortex_case_unstructured
   end type

contains

   subroutine eval_vortex_case(case,nx,ny,rho,ux,uy)
      class(vortex_case_t), intent(in) :: case
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: rho(:,:)
      real(wp), intent(out) :: ux(:,:), uy(:,:)

      real(wp) :: xl, yl, xr, yr, rsqr
      real(wp) :: Rcsqr, vMa_sq
      real(wp), parameter :: half = 1._wp/2._wp
      integer :: x, y

      !$omp parallel default(private) shared(case,nx,ny,rho,ux,uy)

      Rcsqr = case%Rc**2

      ! vortex Mach number
      vMa_sq = case%eps**2 / case%csqr

      !$omp do schedule(static)
      do x = 1, nx
         xl = real(x-1,wp) + 0.5_wp
         do y = 1, ny
            yl = real(y-1,wp) + 0.5_wp

            xr = xl - case%xc
            yr = yl - case%yc

            rsqr = xr*xr + yr*yr

            ! Convected vortex, Eqs. (3) and (4)
            ux(y,x) = case%U0 - case%eps * (yr/case%Rc) * exp(-half*rsqr/Rcsqr)
            uy(y,x) =           case%eps * (xr/case%Rc) * exp(-half*rsqr/Rcsqr)

            ! Barotopic density initialization, Eq. (20)
            rho(y,x) = case%rho0 * exp(-half * vMa_sq * exp(-rsqr/Rcsqr))

            ! First order Taylor expansion, Eq. (21)
            !rho(y,x) = case%rho0 * (1 - half * vMa_sq * exp(-rsqr/Rcsqr))

            ! Second order Taylor expansion, Eq. (22)
            !rho(y,x) = case%rho0 * (1 - half * vMa_sq * exp(-rsqr/Rcsqr) + &
            !   0.125_wp * vMa_sq**2 / gamma * exp(-half*rsqr/Rcsqr))

         end do
      end do
      !$omp end do 

      !$omp end parallel

   end subroutine

   subroutine eval_vortex_case_unstructured(case,n,xy,rho,ux,uy)
      class(vortex_case_t), intent(in) :: case
      integer, intent(in) :: n
      real(wp), intent(in) :: xy(:,:) ! shape [2,n]
      real(wp), intent(out) :: rho(:) ! shape [n]
      real(wp), intent(out) :: ux(:), uy(:) ! shape[n]

      real(wp) :: xr, yr, rsqr
      real(wp) :: Rcsqr, vMa_sq
      real(wp), parameter :: half = 1._wp/2._wp
      integer :: i

      !$omp parallel default(private) shared(case,n,xy,rho,ux,uy)

      Rcsqr = case%Rc**2

      ! vortex Mach number
      vMa_sq = case%eps**2 / case%csqr

      !$omp do schedule(static)
      do i = 1, n

         xr = xy(1,i) - case%xc
         yr = xy(2,i) - case%yc

         rsqr = xr*xr + yr*yr

         ! Convected vortex, Eqs. (3) and (4)
         ux(i) = case%U0 - case%eps * (yr/case%Rc) * exp(-half*rsqr/Rcsqr)
         uy(i) =           case%eps * (xr/case%Rc) * exp(-half*rsqr/Rcsqr)

         ! Barotopic density initialization, Eq. (20)
         rho(i) = case%rho0 * exp(-half * vMa_sq * exp(-rsqr/Rcsqr))

         ! First order Taylor expansion, Eq. (21)
         !rho(i) = case%rho0 * (1 - half * vMa_sq * exp(-rsqr/Rcsqr))

         ! Second order Taylor expansion, Eq. (22)
         !rho(i) = case%rho0 * (1 - half * vMa_sq * exp(-rsqr/Rcsqr) + &
         !   0.125_wp * vMa_sq**2 / gamma * exp(-half*rsqr/Rcsqr))

      end do
      !$omp end do 

      !$omp end parallel

   end subroutine

end module
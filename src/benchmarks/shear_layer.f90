module shear_layer_case

   use precision, only: wp

   implicit none
   private

   public :: shear_layer_t

   type :: shear_layer_t
      real(wp) :: length
      real(wp) :: U0
      real(wp) :: k
      real(wp) :: delta
      real(wp) :: rho0 = 1.0_wp
   contains
      procedure :: eval => eval_shear_layer_unstructured, &
                           eval_shear_layer_grid
   end type

   real(wp), parameter :: pi = 4.0_wp*atan(1.0_wp)

contains

   subroutine eval_shear_layer_grid(case,nx,ny,rho,ux,uy)
      class(shear_layer_t), intent(in) :: case
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: rho(:,:), ux(:,:), uy(:,:)

      real(wp) :: xd, yd
      integer :: x, y

      do x = 1, nx
         ! Scale x-coordinate to region [0,1)
         xd = (real(x-1,wp) + 0.5_wp)/case%length
         do y = 1, ny
            ! Scale y-coordinate to region [0,1)
            yd = (real(y-1,wp) + 0.5_wp)/case%length

            uy(y,x) = case%U0 * case%delta * sin(2*pi*(xd + 0.25_wp))
            if (yd <= 0.5_wp) then
               ux(y,x) = case%U0 * tanh(case%k * (yd - 0.25_wp))
            else
               ux(y,x) = case%U0 * tanh(case%k * (0.75_wp - yd))
            end if

            rho(y,x) = case%rho0

         end do
      end do

   end subroutine

   subroutine eval_shear_layer_unstructured(case,n,xy,rho,ux,uy)
      class(shear_layer_t), intent(in) :: case
      integer, intent(in) :: n
      real(wp), intent(in) :: xy(2,n)
      real(wp), intent(out) :: rho(n), ux(n), uy(n)

      real(wp) :: xd, yd
      integer :: i

      do i = 1, n

         ! Scale coordinate to region [0,1)
         xd = xy(1,i)/case%length
         yd = xy(2,i)/case%length

         uy(i) = case%U0 * case%delta * sin(2*pi*(xd + 0.25_wp))

         if (yd <= 0.5_wp) then
            ux(i) = case%U0 * tanh(case%k * (yd - 0.25_wp))
         else
            ux(i) = case%U0 * tanh(case%k * (0.75_wp - yd))
         end if

         rho(i) = case%rho0

      end do

   end subroutine

end module
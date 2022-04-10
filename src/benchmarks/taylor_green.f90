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

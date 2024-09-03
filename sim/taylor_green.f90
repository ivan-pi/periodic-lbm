module taylor_green
   use kinds, only: wp
   use lbm_primitives, only: fluid_state

   implicit none
   private

   public :: taylor_green_vortex
   public :: tgv_decay_time, tgv_state
   public :: tgv_gnuplot

   public :: pi

   type :: taylor_green_vortex
      real(wp) :: k1, k2, umax, nu
   end type

   real(wp), parameter :: pi = 4._wp*atan(1._wp)

contains

   subroutine tgv_gnuplot(case,filename)
      type(taylor_green_vortex), intent(in) :: case
      character(len=*), intent(in) :: filename
      
      integer :: unit

      open(newunit=unit,file=trim(filename),action="write")

      write(unit,'(A,G0)') "k1 = ", case%k1
      write(unit,'(A,G0)') "k2 = ", case%k2
      write(unit,'(A,G0)') "umax = ", case%umax
      write(unit,'(A,G0)') "nu = ", case%nu
      write(unit,'(A)') "td = 1.0/(nu*(k1**2 + k2**2))"
      write(unit,'(A)') "f(x) = umax*exp(-x/td)"

      close(unit)
   end subroutine

   ! Time constant for vortex decay
   !
   ! Essentially the vortex half-life
   ! 
   pure function tgv_decay_time(case) result(td)
      type(taylor_green_vortex), intent(in) :: case
      real(wp) :: td
      td = decay_time(case%nu,case%k1,case%k2)

   contains

      pure function decay_time(nu,k1,k2)
         real(wp), intent(in) :: nu, k1, k2
         real(wp) :: decay_time
         decay_time = 1.0_wp/(nu*(k1**2 + k2**2))
      end function

   end function

   pure function tgv_state(case,x,y,t) result(state)
      type(taylor_green_vortex), intent(in) :: case
      real(wp), intent(in) :: x, y, t
      type(fluid_state) :: state

      real(wp) :: td

      td = tgv_decay_time(case)

      associate(umax => case%umax, k1 => case%k1, k2 => case%k2)

      state%p = -0.25_wp*(umax**2)*((k2/k1)*cos(2.0_wp*k1*x) + (k1/k2)*cos(2.0_wp*k2*y))*exp(-2.0_wp*t/td)

      state%u = -umax*sqrt(k2/k1)*cos(k1*x)*sin(k2*y)*exp(-t/td)
      state%v =  umax*sqrt(k1/k2)*sin(k1*x)*cos(k2*y)*exp(-t/td)

      ! TODO: double check the shear formulas

      state%sxx = umax*sqrt(k1*k2)*sin(k1*x)*sin(k2*y)*exp(-t/td)
      state%sxy = 0.5_wp*umax*(sqrt(k1**3/k2) - sqrt(k2**3/k1))*cos(k1*x)*cos(k2*y)*exp(-t/td)
      state%syy = -state%sxx

      end associate

   end function

end module

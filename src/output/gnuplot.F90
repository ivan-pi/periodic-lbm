module output_gnuplot

   use precision, only: wp
   implicit none
   private

   public :: output_gnuplot_grid

#if PRECISION_SP
   !> Format string for single precision real numbers
   character(*), parameter :: FMT_REAL = '(ES15.8E2)'
#else
   !> Format string for double precision real numbers
   character(*), parameter :: FMT_REAL = '(ES24.16E3)'
#endif


contains

   subroutine output_gnuplot_grid(filename,nx,ny,rho,ux,uy)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: rho(ny,nx), ux(ny,nx), ux(ny,nx)

      integer :: unit, x, y
      real(wp) :: xx, yy

      open(newunit=unit,file=filename)

      do x = 1, grid%nx
         xx = (x - 1) + 0.5_wp
         do y = 1, grid%ny
            yy = (y - 1) + 0.5_wp
            write(unit,'(5('//FMT_REAL//',:,1X))') xx, yy, grid%rho(y,x), grid%ux(y,x), grid%uy(y,x)
         end do
         write(unit,*)
      end do

      close(unit)
   end subroutine
   
end module
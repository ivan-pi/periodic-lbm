module output_npy

   use precision, only: wp
   use stdlib_io_npy, only: save_npy

   implicit none
   private

   public :: output_fluid_npy

contains

   subroutine output_fluid_npy(filename,nx,ny,mf)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx, ny
      real(wp), intent(in), target :: mf(ny,nx,3)

      integer :: iostat
      character(len=:), allocatable :: iomsg

      call save_npy(filename, mf, iostat, iomsg)
      if (iostat /= 0) then
         print *, iomsg
      end if

   end subroutine

end module

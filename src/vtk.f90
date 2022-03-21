module vtk

   implicit none
   private

   public :: output_vtk_grid_ascii
   public :: output_vtk_grid_binary

contains

   subroutine output_vtk_grid_ascii(filename,nx,ny,rho,ux,uy)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx, ny
      real(wp), intent(in), dimension(ny,nx) :: rho, ux, uy

      integer :: unit_, i, j 

      open(newunit=unit_,file=filename)


      write(unit_, "(A)") "# vtk DataFile Version 3.0"
      write(unit_, "(A)") "fluid"
      write(unit_, "(A)") "ASCII"
      write(unit_, "(A)") "DATASET RECTILINEAR_GRID"
      write(unit_, "(A,I5,I5,I5)") "DIMENSIONS ", nx, ny, 1
      write(unit_, "(A,I5,A)") "X_COORDINATES ", nx, " float"
      do i = 1, nx
         write(unit_, "(ES12.5)") real(i-1,wp)
      end do
      write(unit_, "(A,I5,A)") "Y_COORDINATES ", ny, " float"
      do j = 1, ny
         write(unit_, "(ES12.5)") real(j-1,wp)
      end do
      write(unit_, "(A,I5,A)") "Z_COORDINATES ", 1, " float"
      write(unit_, "(ES12.5)") 0.0_wp
      write(unit_, "(A,I5)") "POINT_DATA ", nx*ny
      write(unit_, "(A)") "SCALARS density float 1"
      write(unit_, "(A)") "LOOKUP_TABLE default"
      do j = 1, ny
         do i = 1, nx
            write(unit_, "(ES12.5)") rho(i,j)
         end do
      end do
      write(unit_, *) "VECTORS velocity float"
      do j = 1, ny
         do i = 1, nx
            write(unit_, "(3ES12.5)") ux(i,j), uy(i,j), 0.0_wp
         end do
      end do

      close(unit_)   

   end subroutine output_vtk_grid_ascii

   subroutine output_vtk_grid_binary(filename,nx,ny,rho,ux,uy)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: rho(ny,nx)
      real(wp), intent(in), target :: ux(ny,nx), uy(ny,nx)

      real(wp), pointer :: uxf(:) => null()
      real(wp), pointer :: uyf(:) => null()
      character(len=8) :: nxstr, nystr, nzstr, nstr
      character(len=1), parameter :: nl = new_line('a')
      integer :: vtkunit, i

      associate(n => nx*ny)

      ! Flatten 2D arrays using pointer bound remapping
      uxf(1:n) => ux
      uyf(1:n) => uy

      nxstr = i2str(nx)
      nystr = i2str(ny)
      nzstr = i2str(1)

      nstr = i2str(n)


      open(newunit=vtkunit,file=filename,access="stream",convert="BIG_ENDIAN")


      write(vtkunit) '# vtk DataFile Version 3.0' // nl
      write(vtkunit) 'fluid' // nl
      write(vtkunit) 'BINARY' // nl
      write(vtkunit) 'DATASET RECTILINEAR_GRID' // nl
      write(vtkunit) 'DIMENSIONS ' // nxstr // nstr // nzstr // nl
      write(vtkunit) 'X_COORDINATES ' // nxstr // ' DOUBLE' // nl
      write(vtkunit) (real(i-1,wp), i = 1, nx), nl
      write(vtkunit) 'Y_COORDINATES ' // nystr //' DOUBLE' // nl
      write(vtkunit) (real(i-1,wp), i = 1, ny)
      write(vtkunit) 'Z_COORDINATES ' // nzstr // ' DOUBLE' // nl
      write(vtkunit) 0.0_wp, nl
      write(vtkunit) "POINT_DATA", nstr, nl
      write(vtkunit) "SCALARS density DOUBLE 1"//nl
      write(vtkunit) "LOOKUP_TABLE default", nl
      write(vtkunit) rho, nl
      write(vtkunit) "VECTORS velocity DOUBLE 3"//nl
      write(vtkunit) (uxf(i), uyf(i), 0, i = 1, n),nl

      close(vtkunit)        
     
      end associate

   contains

      pure function i2str(i) result(str)
         integer, intent(in) :: i
         character(len=8) :: str
         write(str,'(i8)') i     
      end function i2str

   end subroutine output_vtk_grid_binary



end module
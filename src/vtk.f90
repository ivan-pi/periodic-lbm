module vtk

   implicit none
   private

   public :: output_vtk_grid_ascii
   public :: output_vtk_grid_binary

   public :: output_vtk_structuredPoints


   integer, parameter :: wp = kind(1.0d0)
   character(*), parameter :: FMT_REAL = '(es15.8e2)'

   interface
! void output_vtk_polydata(
!     const char*, 
!     int, 
!     const double*, 
!     const double*, 
!     const double*, 
!     const double*)
      subroutine c_output_vtk_polydata(filename,n,p,rho,ux,uy) bind(c)
         use iso_c_binding, only: c_char, c_double, c_int
         implicit none
         character(kind=c_char), intent(in) :: filename(*)
         integer(c_int), intent(in), value :: n
         real(c_double), intent(in) :: p(2,n)
         real(c_double), intent(in) :: rho(n)
         real(c_double), intent(in) :: ux(n)
         real(c_double), intent(in) :: uy(n)
      end subroutine
   end interface

contains

   subroutine output_vtk_grid_ascii(filename,nx,ny,rho,ux,uy)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx, ny
      real(wp), intent(in), dimension(ny,nx) :: rho, ux, uy

      integer :: unit_, i, j 

      open(newunit=unit_,file=filename)


      write(unit_, '(A)') "# vtk DataFile Version 3.0"
      write(unit_, '(A)') "fluid"
      write(unit_, '(A)') "ASCII"
      write(unit_, '(A)') "DATASET RECTILINEAR_GRID"
      write(unit_, '(A,3(I0,:,1X))') "DIMENSIONS ", nx+1, ny+1, 2
      write(unit_, '(A,I0,A)') "X_COORDINATES ", nx+1, " float"
      do i = 1, nx+1
         write(unit_, '(ES12.5)') real(i-1,wp)
      end do
      write(unit_, '(A,I0,A)') "Y_COORDINATES ", ny+1, " float"
      do j = 1, ny+1
         write(unit_, '(ES12.5)') real(j-1,wp)
      end do
      write(unit_, '(A)') "Z_COORDINATES 2 float"
      write(unit_, '(ES12.5)') 0.0_wp
      write(unit_, '(ES12.5)') 1.0_wp
      write(unit_, '(A,I5)') "CELL_DATA ", nx*ny
      write(unit_, '(A)') "SCALARS density float 1"
      write(unit_, '(A)') "LOOKUP_TABLE default"
      do i = 1, nx
         do j = 1, ny
            write(unit_, '(ES12.5)') rho(j,i)
         end do
      end do
      write(unit_, *) "VECTORS velocity float"
      do i = 1, nx
         do j = 1, ny
            write(unit_, '(3ES12.5)') ux(j,i), uy(j,i), 0.0_wp
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




   subroutine output_vtk_structuredPoints(filename,nx,ny,rho,ux,uy)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx, ny
      real(wp), intent(in), dimension(ny,nx) :: rho, ux, uy

      integer :: unit_, i, j 

      open(newunit=unit_,file=filename)


      write(unit_, '(A)') "# vtk DataFile Version 3.0"
      write(unit_, '(A)') "fluid"
      write(unit_, '(A)') "ASCII"
      write(unit_, '(A)') "DATASET STRUCTURED_POINTS"
      write(unit_, '(A,3(I0,1X))') "DIMENSIONS ", nx+1, ny+1, 2
      write(unit_, '(A,3(ES15.8E2,1X))') "ORIGIN ", 0.0_wp, 0.0_wp, 0.0_wp
      write(unit_, '(A,3(ES15.8E2,1X))') "SPACING ", 1.0_wp, 1.0_wp, 1.0_wp
      
      !
      ! Dataset attributes
      !
      write(unit_, *)
      write(unit_, '(A,I0)') "CELL_DATA ", nx*ny

      write(unit_, '(A)') "SCALARS Density float 1"
      write(unit_, '(A)') "LOOKUP_TABLE default"
      do i = 1, nx
         do j = 1, ny
            write(unit_, '(ES15.8E2)') rho(j,i)
         end do
      end do

      write(unit_, *)
      write(unit_, '(A)') "VECTORS Velocity float"
      do i = 1, nx
         do j = 1, ny
            write(unit_, '(3(ES15.8E2,1X))') ux(j,i), uy(j,i), 0.0_wp
         end do
      end do

      close(unit_)   

   end subroutine output_vtk_structuredPoints


end module
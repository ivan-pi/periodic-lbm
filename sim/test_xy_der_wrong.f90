! test_xy_der.f90 -- 
!     Test for effective bandwidth of 2-d mixed derivative FD stencil
!
!     Uses OpenACC, stdpar, OpenMP, and CUDA Fortran; in other words
!     don't forget the flags `-cuda -acc=gpu -mp=gpu -stdpar=gpu`
!
! Ivan Pribec, September 3, 2024

module xy_der_mod

implicit none
public

integer, parameter :: wp = kind(1.0d0)

contains

   subroutine set_known_field(nx,ny,f,f_xy)
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: f(nx,ny)
      real(wp), intent(out) :: f_xy(nx,ny)

      real(wp), parameter :: pi = 4*atan(1.0_wp)
      real(wp) :: x, y, dx, dy
      integer :: i, j

      do j = 1, ny
         do i = 1, nx

            x = (i-1) + 0.5_wp
            y = (j-1) + 0.5_wp

            ! known function
            f(i,j) = sin(2*pi*(x/nx))*cos(2*pi*(y/ny))

            ! second partial derivative w.r.t x- and y-coordinates
            f_xy(i,j) = (-4*pi**2/(real(nx,wp)*real(ny,wp))) &
               * cos(2*pi*(x/ny))*sin(2*pi*(y/ny))

         end do
      end do

   end subroutine

   ! do concurrent version
   subroutine cross_derivative_dc(nx,ny,f,f_xy)
      integer, intent(in), value :: nx, ny
      real(wp), intent(in) :: f(0:nx+1,0:ny+1) ! includes HALO layer
      real(wp), intent(out) :: f_xy(nx,ny)

      integer :: i, j

      do concurrent(i=1:nx,j=1:ny)
         f_xy(i,j) = 0.25_wp*(f(i+1,j+1) - f(i-1,j+1) + f(i-1,j-1) - f(i+1,j-1))
      end do

   end subroutine

   ! OpenACC version
   subroutine cross_derivative_acc(nx,ny,f,f_xy)
      integer, intent(in), value :: nx, ny
      real(wp), intent(in) :: f(0:nx+1,0:ny+1) ! includes HALO layer
      real(wp), intent(out) :: f_xy(nx,ny)

      integer :: i, j

      !$acc parallel loop collapse(2) copyin(f) copyout(f_xy)
      do j = 1, ny
         do i = 1, nx
            f_xy(i,j) = 0.25_wp*(f(i+1,j+1) - f(i-1,j+1) + f(i-1,j-1) - f(i+1,j-1))
         end do
      end do

   end subroutine

   ! OpenMP version
   subroutine cross_derivative_omp(nx,ny,f,f_xy)
      integer, intent(in), value :: nx, ny
      real(wp), intent(in) :: f(0:nx+1,0:ny+1) ! includes HALO layer
      real(wp), intent(out) :: f_xy(nx,ny)

      integer :: i, j

      !$omp target teams distribute parallel do collapse(2) map(to: f) map(from: f_xy)
      do j = 1, ny
         do i = 1, nx
            f_xy(i,j) = 0.25_wp*(f(i+1,j+1) - f(i-1,j+1) + f(i-1,j-1) - f(i+1,j-1))
         end do
      end do

   end subroutine

   ! CUDA Fortran kernel routine
   attributes(global) subroutine cross(nx,ny,f,f_xy)
      integer, intent(in), value :: nx, ny
      real(wp), intent(in) :: f(0:nx+1,0:ny+1) ! includes HALO layer
      real(wp), intent(out) :: f_xy(nx,ny)

      integer :: i, j

      i = (blockIdx%x - 1)*blockDim%x + threadIdx%x
      j = (blockIdx%y - 1)*blockDim%y + threadIdx%y

      if (i <= nx .and. j <= ny) then
         f_xy(i,j) = 0.25_wp*(f(i+1,j+1) - f(i-1,j+1) + f(i-1,j-1) - f(i+1,j-1))
      end if

   end subroutine

   ! CUDA Fortran version
   subroutine cross_derivative_cuf(nx,ny,f,f_xy)
      
      use cudafor, only: dim3

      integer, intent(in), value :: nx, ny
      real(wp), intent(in), device :: f(0:nx+1,0:ny+1) ! includes HALO layer
      real(wp), intent(out), device :: f_xy(nx,ny)

      type(dim3) :: grid, tBlock
      
      tBlock = dim3(64,8,1)
      grid = dim3((nx + tBlock%x - 1)/tBlock%x, &
                  (nx + tBlock%y - 1)/tBlock%y, 1)

      call cross<<<grid,tBlock>>>(nx,ny,f,f_xy)

   end subroutine

end module


program xy_der_test

use xy_der_mod
implicit none

integer :: nx, ny, nargs

real(wp), allocatable :: f(:,:), f_xy(:,:), approx_xy(:,:)
real(wp), allocatable, device :: fd(:,:), fd_xy(:,:)
integer :: i, j

nargs = command_argument_count()
if (nargs /= 2) then
   call print_usage
   stop
end if

call parse_args(nx,ny)
call print_gpu_info

allocate(f(nx,ny), f_xy(nx,ny))
allocate(fd(0:nx+1,0:ny+1))

allocate(fd_xy(nx,ny))
allocate(approx_xy(nx,ny))

! Set the field to a known periodic function
call set_known_field(nx,ny,f,f_xy)

! Copy Host to Device
fd(1:nx,1:ny) = f

! Fill HALO layers (SOUTH and NORTH)
do concurrent(i=1:nx)
   fd(i,0)    = fd(i,ny)
   fd(i,ny+1) = fd(i,1)
end do

! Fill HALO layers (EAST and WEST)
do concurrent(j=0:ny+1)
   fd(0,j) = fd(ny,j)
   fd(nx+1,j) = fd(1,j)
end do

call measure(cross_derivative_dc, "Do-Concurrent")
call measure(cross_derivative_acc, "OpenACC-parallel-loop")
call measure(cross_derivative_acc, "OpenMP-target-loop")
call measure(cross_derivative_cuf, "CUDA-Fortran-naive")

contains

   subroutine print_usage
      write(*,'(A)') "Usage: test_xy_der <nx> <ny>"
   end subroutine

   subroutine parse_args(nx,ny)
      integer, intent(out) :: nx, ny
      character(len=64) :: str
      call get_command_argument(1,str)
      read(str,*) nx
      call get_command_argument(2,str)
      read(str,*) ny
   end subroutine

   subroutine print_gpu_info
      use cudafor
      type(cudaDeviceProp) :: prop
      integer :: istat

      istat = cudaGetDeviceProperties(prop,0)

      print '("Device name: ", A)', trim(prop%name)
      print '("  Compute Capability: ",I0,".",I0)', prop%major, prop%minor
      print '("  Memory Clock Rate (KHz): ", I0)', prop%memoryClockRate
      print '("  Memory Bus Width (bits): ", I0)', prop%memoryBusWidth
      print '("  Peak Memory Bandwith (GB/s): ", f9.2)', &
         2.0 * prop%memoryClockRate * (prop%memoryBusWidth / 8) * 1.0e-6

      print '()'
   end subroutine

   ! Measurement helper procedure
   subroutine measure(kernel,name)
      use cudafor
      interface
         subroutine kernel(nx,ny,fd,fd_xy)
            import wp
            integer, intent(in), value :: nx, ny
            real(wp), intent(in), device :: fd(0:nx+1,0:ny+1)
            real(wp), intent(out), device :: fd_xy(nx,ny)
         end subroutine
      end interface
      character(len=*), intent(in) :: name
      
      integer :: istat
      type(cudaEvent) :: startEvent, stopEvent
      real(kind(1.0e0)) :: time, elapsed

      integer :: niter, k


      istat = cudaEventCreate(startEvent)
      istat = cudaEventCreate(stopEvent)
      niter = 1

      outer: do
         istat = cudaEventRecord(startEvent,0)
         do k = 1, niter
            call kernel(nx,ny,fd,fd_xy)
         end do
         istat = cudaEventRecord(stopEvent,0)
         istat = cudaEventSynchronize(stopEvent)
         istat = cudaEventElapsedTime(time,startEvent,stopEvent)

         ! Aim for at least 0.2 s total runtime
         if (time > 200.0) exit outer
         niter = 2*niter
      end do outer

      istat = cudaEventDestroy(startEvent)
      istat = cudaEventDestroy(stopEvent)

      ! Copy results to host 
      approx_xy = fd_xy

      elapsed = time/niter

      write(*,*) trim(name)
      write(*,*) "-----------------------------------------------------v"
      write(*,*) "Avg. time (ms): ", elapsed
      write(*,*) "Max. absolute error: ", maxval(abs(f_xy - approx_xy))
      write(*,*) "Sum of absolute residuals: ", sum(abs(f_xy - approx_xy))
      write(*,*) "Sum of squared residuals: ", sum((f_xy - approx_xy)**2)
      write(*,*) "Effective Bandwidth (GB/s): ", bw(nx,ny,elapsed/1000.0)
      write(*,*)

   end subroutine

   ! Effective bandwidth helper function
   pure function bw(nx,ny,time)
      integer, intent(in) :: nx, ny
      real(kind(1.0e0)), intent(in) :: time ! Time in seconds
      real(kind(1.0e0)) :: rb, wb, bw
      rb = real(nx+1,wp)*real(ny+1,wp)*8
      wb = real(nx,wp)*real(ny,wp)*8
      bw = ((rb + wb) * 1.0e-9) / time
   end function

end program
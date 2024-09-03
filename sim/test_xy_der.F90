! test_xy_der.F90 -- 
!     Test for effective bandwidth of 2-d mixed derivative FD stencil
!
! Ivan Pribec, September 3, 2024
!
module xy_der_mod

implicit none
public

integer, parameter :: wp = kind(1.0d0)

contains

   subroutine output(nx,ny,f1,f2)
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: f1(nx,ny), f2(nx,ny)

      real(wp) :: x, y
      integer :: i, j, wunit
      character(len=60) :: fname

      fname = "output.dat"

      open(newunit=wunit,file=trim(fname),status="unknown",action="write")

      do j = 1, ny
         do i = 1, nx

            x = (i-1) + 0.5_wp
            y = (j-1) + 0.5_wp

            write(wunit,'(4(ES24.17,:,1X))') x, y, f1(i,j), f2(i,j)

         end do
      end do

      close(wunit)

   end subroutine

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
   subroutine cross_derivative_dc(nx,ny,ldx,f,f_xy)
      integer, intent(in), value :: nx, ny,ldx
      real(wp), intent(in) :: f(0:nx+1,0:ny+1) ! includes HALO layer
      real(wp), intent(out) :: f_xy(ldx,ny)

      integer :: i, j

      do concurrent(i=1:nx,j=1:ny)
         f_xy(i,j) = 0.25_wp*(f(i+1,j+1) - f(i-1,j+1) + f(i-1,j-1) - f(i+1,j-1))
      end do

   end subroutine

   ! OpenACC version
   subroutine cross_derivative_acc(nx,ny,ldx,f,f_xy)
      integer, intent(in), value :: nx, ny, ldx
      real(wp), intent(in) :: f(0:nx+1,0:ny+1) ! includes HALO layer
      real(wp), intent(out) :: f_xy(ldx,ny)

      integer :: i, j

      !$acc parallel loop collapse(2)
      do j = 1, ny
         do i = 1, nx
            f_xy(i,j) = 0.25_wp*(f(i+1,j+1) - f(i-1,j+1) + f(i-1,j-1) - f(i+1,j-1))
         end do
      end do

   end subroutine

   ! OpenMP version
   subroutine cross_derivative_omp(nx,ny,ldx,f,f_xy)
      integer, intent(in), value :: nx, ny, ldx
      real(wp), intent(in) :: f(0:nx+1,0:ny+1) ! includes HALO layer
      real(wp), intent(out) :: f_xy(ldx,ny)

      integer :: i, j

      !$omp target teams distribute parallel do collapse(2)
      do j = 1, ny
         do i = 1, nx
            f_xy(i,j) = 0.25_wp*(f(i+1,j+1) - f(i-1,j+1) + f(i-1,j-1) - f(i+1,j-1))
         end do
      end do

   end subroutine

   ! CUDA Fortran kernel routine
   attributes(global) subroutine cross(nx,ny,ldx,f,f_xy)
      integer, intent(in), value :: nx, ny, ldx
      real(wp), intent(in) :: f(0:nx+1,0:ny+1) ! includes HALO layer
      real(wp), intent(out) :: f_xy(ldx,ny)

      integer :: i, j

      i = (blockIdx%x - 1)*blockDim%x + threadIdx%x
      j = (blockIdx%y - 1)*blockDim%y + threadIdx%y

      if (i <= nx .and. j <= ny) then
         f_xy(i,j) = 0.25_wp*(f(i+1,j+1) - f(i-1,j+1) + f(i-1,j-1) - f(i+1,j-1))
      end if

   end subroutine

   ! CUDA Fortran version
   subroutine cross_derivative_cuf(nx,ny,ldx,f,f_xy)
      
      use cudafor, only: dim3

      integer, intent(in), value :: nx, ny, ldx
      real(wp), intent(in), device :: f(0:nx+1,0:ny+1) ! includes HALO layer
      real(wp), intent(out), device :: f_xy(ldx,ny)

      type(dim3) :: grid, tBlock
      
      tBlock = dim3(64,8,1)
      grid = dim3((nx + tBlock%x - 1)/tBlock%x, &
                  (nx + tBlock%y - 1)/tBlock%y, 1)

      call cross<<<grid,tBlock>>>(nx,ny,ldx,f,f_xy)

   end subroutine

#define BLOCK_X 32
#define BLOCK_Y 8

   ! CUDA Fortran kernel routine
   attributes(global) subroutine cross_shm(nx,ny,f,f_xy)
      integer, intent(in), value :: nx, ny
      real(wp), intent(in) :: f(0:nx+1,0:ny+1) ! includes HALO layer
      real(wp), intent(out) :: f_xy(nx,ny)

      integer, parameter :: radius = 1
      real(wp), shared :: fs(1-radius:BLOCK_X+radius,1-radius:BLOCK_Y+radius)

      integer :: i, j, is, js

      i = (blockIdx%x - 1)*blockDim%x + threadIdx%x
      j = (blockIdx%y - 1)*blockDim%y + threadIdx%y

      is = threadIdx%x
      js = threadIdx%y

      ! The shared memory block is bigger than the block dimensions,
      ! meaning some threads need to do the extra work of filling the halo

      if (i <= nx+radius .and. j <= ny+radius) fs(is,js) = f(i,j)

      if (is <= radius .and. j <= ny+radius .and. i <= nx) then
         fs(is - radius,js) = f(i - radius,j)
         fs(is + BLOCK_X,js) = f(i + BLOCK_X,j)
      end if
      if (js <= radius .and. i <= nx+radius .and. j <= ny) then
         fs(is,js - radius) = f(i,j - radius)
         fs(is,js + BLOCK_Y) = f(i,j + BLOCK_Y)
      end if
      if (i <= nx .and. j <= ny .and. is <= radius .and. js <= radius) then
         fs(is - radius,js-radius) = f(i - radius,j - radius)
         fs(is - radius,js+BLOCK_Y) = f(i-radius,j+BLOCK_Y)
         fs(is + BLOCK_X,js+BLOCK_Y) = f(i+BLOCK_X,j+BLOCK_Y)
         fs(is + BLOCK_X,js-radius) = f(i+BLOCK_X,j - radius)
      end if
      call syncthreads

      if (i <= nx .and. j <= ny) then
         f_xy(i,j) = 0.25_wp*(fs(is+1,js+1) - fs(is-1,js+1) + fs(is-1,js-1) - fs(is+1,js-1))
      end if

   end subroutine

   ! CUDA Fortran version
   subroutine cross_derivative_cufshm(nx,ny,f,f_xy)
      
      use cudafor, only: dim3

      integer, intent(in), value :: nx, ny
      real(wp), intent(in), managed :: f(0:nx+1,0:ny+1) ! includes HALO layer
      real(wp), intent(out), managed :: f_xy(nx,ny)

      type(dim3) :: grid, tBlock
      
      tBlock = dim3(BLOCK_X,BLOCK_Y,1)
      grid = dim3((nx + tBlock%x - 1)/tBlock%x, &
                  (nx + tBlock%y - 1)/tBlock%y, 1)

      call cross_shm<<<grid,tBlock>>>(nx,ny,f,f_xy)

   end subroutine

end module


program xy_der_test

use xy_der_mod
implicit none

integer :: ldx, nx, ny, nargs

real(wp), allocatable :: f(:,:), f_xy(:,:), approx_xy(:,:)
real(wp), allocatable :: fd(:,:), fd_xy(:,:)
integer :: i, j

nargs = command_argument_count()
if (nargs /= 2) then
   call print_usage
   stop
end if

call parse_args(nx,ny)
call print_gpu_info

! Pad to a multiple of the warp size
ldx = ((nx + 64 - 1) / 64) * 64

print *, "Problem dimensions"
print *, "  ldx = ", ldx
print *, "  nx  = ", nx
print *, "  ny  = ", ny
print *, ""

allocate(f(nx,ny), f_xy(nx,ny))
allocate(fd(0:nx+1,0:ny+1))

allocate(fd_xy(ldx,ny))
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

!call measure(cross_derivative_cufshm, "CUDA-Fortran-shared-mem")
!call output(nx,ny,f_xy,approx_xy)

!deallocate(fd,fd_xy)
!deallocate(f,f_xy,approx_xy)

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
         subroutine kernel(nx,ny,ldx,fd,fd_xy)
            import wp
            integer, intent(in), value :: nx, ny, ldx
            real(wp), intent(in), device :: fd(0:nx+1,0:ny+1)
            real(wp), intent(out), device :: fd_xy(ldx,ny)
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
            call kernel(nx,ny,ldx,fd,fd_xy)
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
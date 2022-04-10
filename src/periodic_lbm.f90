module periodic_lbm

   use precision, only: wp
   use fvm_bardow, only: lattice_grid, cx, cy

   implicit none
   private

   public :: perform_lbm_step

   public :: lbm_stream

contains

   subroutine perform_lbm_step(grid)
      type(lattice_grid), intent(inout) :: grid


      call grid%streaming() ! write from iold to inew 
      call grid%collision() ! update inew in place

      block
         integer :: itmp
         itmp = grid%iold
         grid%iold = grid%inew
         grid%inew = itmp
      end block

   end subroutine


   subroutine lbm_stream(grid)
      class(lattice_grid), intent(inout) :: grid

      integer :: ld 

      ld = size(grid%f,1)

      call lbm_stream_kernel(grid%nx,grid%ny,ld, &
         fsrc = grid%f(:,:,:,grid%iold), &
         fdst = grid%f(:,:,:,grid%inew))

   contains

      subroutine lbm_stream_kernel(nx,ny,ld,fsrc,fdst)
         integer, intent(in) :: nx, ny, ld
         real(wp), intent(in) :: fsrc(ld,nx,0:8)
         real(wp), intent(out) :: fdst(ld,nx,0:8)

         integer :: x, y
         integer :: xp1, xm1, yp1, ym1


      !$omp parallel default(private) shared(nx,ny,ld,fsrc,fdst)
      
      !$omp workshare
      fdst(1:ny,:,0) = fsrc(1:ny,:,0)
      !$omp end workshare

      !$omp do schedule(static)
      do x = 1, nx

         xp1 = mod(x, nx) + 1
         xm1 = mod(nx + x - 2, nx) + 1

         do y = 1, ny

            fdst(y,x,1) = fsrc(y,xm1,1)
            fdst(y,x,3) = fsrc(y,xp1,3)
         end do
      end do
      !$omp end do

      !$omp do schedule(static)
      do x = 1, nx

         do y = 1, ny

            yp1 = mod(y, ny) + 1
            ym1 = mod(ny + y - 2, ny) + 1

            fdst(y,x,2) = fsrc(ym1,x,2)
            fdst(y,x,4) = fsrc(yp1,x,4)

         end do
      end do
      !$omp end do

      !$omp do schedule(static)
      do x = 1, nx

         xp1 = mod(x, nx) + 1
         xm1 = mod(nx + x - 2, nx) + 1

         do y = 1, ny

            yp1 = mod(y, ny) + 1
            ym1 = mod(ny + y - 2, ny) + 1

            fdst(y,x,5) = fsrc(ym1,xm1,5)
            fdst(y,x,7) = fsrc(yp1,xp1,7)

         end do
      end do
      !$omp end do

      !$omp do schedule(static)
      do x = 1, nx

         xp1 = mod(x, nx) + 1
         xm1 = mod(nx + x - 2, nx) + 1

         do y = 1, ny

            yp1 = mod(y, ny) + 1
            ym1 = mod(ny + y - 2, ny) + 1

            fdst(y,x,6) = fsrc(ym1,xp1,6)
            fdst(y,x,8) = fsrc(yp1,xm1,8)

         end do
      end do
      !$omp end do

      !$omp end parallel

      end subroutine

   end subroutine

end module
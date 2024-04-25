module sparse_streaming

   use, intrinsic :: iso_c_binding
   
   use precision, only: wp
   
   implicit none
   private


   type :: lattice_cloud

      integer :: n, q
      integer :: nnz

      ! Particle distribution function
      real(wp), allocatable :: f(:,:,:)

      ! Macroscopic fields
      ! real(wp), allocatable :: rho(:), ux(:), uy(:)

      ! Sparse matrix storage
      integer, allocatable :: ia(:), ja(:)
      real(wp), allocatable :: lambda(:,:)

      ! MKL handles
      type(sparse_matrix_t), allocatable :: mkl_(:)
      type(matrix_descr) :: descriptor

      ! RSB handles
      type(c_ptr), allocatable :: rsb_(:)

      ! Sparse BLAS handles
      integer, allocatable :: sb_(:)

      ! Indexing
      integer, allocatable :: perm(:)

      ! User-definable procedures
      procedure(collision_interface), pointer :: collide => null()
      procedure(streaming_interface), pointer :: stream  => null()
   
      procedure(initialize_interface), pointer :: initialize => null()

   contains

      procedure :: set_properties
      procedure :: set_relaxation_times
      
   
   end type

contains

   !> Semi-Lagrangian streaming kernel - D2Q9 
   !
   !  Matrix handle is a type(c_ptr)
   !
   subroutine stream_sl_rsb(n,mat,fsrc,fdst)


      use rsb, only: rsb_transposition_n, rsb_idx_kind, &
         rsb_spmv, rsb_error_no_error

      integer, intent(in) :: n
      type(sparse_mat_set_t), intent(in) :: smat
      real(wp), intent(in) :: fsrc(n,0:8)
      real(wp), intent(out) :: fdst(n,0:8)

      integer(c_int), parameter :: op = rsb_transposition_n
      integer(rsb_idx_kind), parameter :: incx = 1, incy = 1

      real(c_double), target :: alpha, beta
      type(c_ptr), target :: mtxp = c_null_ptr

      integer :: q

      alpha = 1.0_c_double
      beta  = 0.0_c_double

      do q = 1, 8

         mtxp = handle(q)
         istat = rsb_spmv(op,c_loc(alpha),mtxp,c_loc(x),incx, &
            c_loc(beta),c_loc(y),incy)

         if (istat /= rsb_error_no_error) then
            ! handle error
         end if

      end do

   end subroutine

   !> Semi-Lagrangian streaming kernel - D2Q9 
   !
   !  Matrix handle is an integer of default kind
   !
   subroutine stream_sl_spblas(n,fsrc,fdst)
      
      use blas_sparse
      use rsb, only: idx_kind => rsb_blas_idx_kind, &
                     ist_kind => rsb_blas_ist_kind
      
      integer, intent(in) :: n
      real(wp), intent(in) :: fsrc(n,0:8)
      real(wp), intent(in) :: fdst(n,0:8)

      integer, parameter :: op = 
      real(wp), parameter :: alpha = 1.0_wp

      integer(idx_kind), parameter :: incx = 1, incy = 1
      integer(ist_kind) :: istat

      integer :: q

      do q = 1, 8
         
         call usmv(op,alpha,handle(q),fsrc(:,q),incx,fdst(:,q),incy,istat)
         
         if (istat /= 0) then
            ! handle error
         end if
     
      end do

   end subroutine


   !> Semi-Lagrangian streaming kernel - D2Q9 
   !
   !  Matrix handle is a type(sparse_matrix_t)
   !  Additionaly a descriptor of type(matrix_descr) is needed
   !  Both are defined in the module mkl_spblas
   !
   subroutine stream_sl_mkl(n,fsrc,fdst)

      use mkl_spblas, only: sparse_operation_non_transpose, &
                            mkl_sparse_d_mv, &
                            sparse_status_success

      integer, intent(in) :: n
      real(wp), intent(in) :: fsrc(n,0:8)
      real(wp), intent(out) :: fdst(n,0:8)

      integer(c_int), parameter :: op = sparse_operation_non_transpose
      real(wp), parameter :: alpha = 1.0_wp, beta = 0.0_wp
      integer(c_int) :: istat
      
      integer :: q

      do q = 1, 8
         
         istat = mkl_sparse_d_mv(op,alpha,handle(q),descr, &
            fsrc,beta,fs(:,q))
      
         if (istat /= sparse_status_success) then
            ! handle error
         end if
      
      end do

   end subroutine


end module
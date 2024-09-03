!> Module implementing standard LBM method
module sim_slbm_class
   
   use, intrinsic :: iso_fortran_env, only: error_unit
   use, intrinsic :: iso_c_binding, only: &
      c_int, c_double, c_ptr, c_loc, c_f_pointer, c_associated
   
   use sim_class, only: sim, dp
   use lbm_primitives, only: &
      lbm_grid, lbm_eqinit, lbm_macros, lbm_collide_and_stream_fused, &
      lbm_periodic_bc_push

   implicit none
   private

   ! Fortran interface
   public :: sim_slbm

   ! Plugin routines
   public :: c_slbm_init
   public :: c_slbm_step
   public :: c_slbm_vars
   public :: c_slbm_free

   type, extends(sim) :: sim_slbm
#if defined(__GFORTRAN__) || defined(__NVCOMPILER)
      type(lbm_grid) :: grid
#else
      type(lbm_grid(:,:)), allocatable :: grid
#endif
      integer :: flip = 0
      real(dp) :: dt = 1.0
   contains
      procedure :: init => slbm_init
      procedure :: step => slbm_step
      procedure :: vars => slbm_vars
   end type

contains

   subroutine slbm_init(this,nx,ny,dt,rho,u,sigma)
      class(sim_slbm), intent(out) :: this
      integer, intent(in) :: nx, ny
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: rho(nx,ny), u(nx,ny,2)    ! rho, (u, v)
      real(dp), intent(in), optional :: sigma(nx,ny,3)  ! (sxx, sxy, syy)

      if (dt /= 1.0_dp) then
         write(error_unit,*) "Standard LBM only supports dt = 1.0!"
         error stop "Failed."
      end if

#if defined(__GFORTRAN__) || defined(__NVCOMPILER)
      this%grid%nx = nx
      this%grid%ny = ny
      allocate(this%grid%f1(0:nx+1,0:ny+1,0:8))
      allocate(this%grid%f2(0:nx+1,0:ny+1,0:8))
#else
      ! PDTs for the win!
      allocate(lbm_grid(nx,ny) :: this%grid)
#endif
      this%flip = 0

      ! Initialize the PDFs using the provided fields
      if (present(sigma)) then
         ! non-equilibrium initialization
!         call lbm_neqinit(nx,ny,this%grid%f1, &
!            rho,u(:,:,1),u(:,:,2), &
!            sigma(:,:,1),sigma(:,:,2),sigma(:,:,3))
         error stop "NotImplementedError"
      else
         ! equilibrium initialization
         call lbm_eqinit(nx,ny,1,this%grid%f1,rho,u(:,:,1),u(:,:,2))
      end if

      this%grid%f2 = this%grid%f1

   end subroutine

   subroutine slbm_step(this,omega)
      class(sim_slbm), intent(inout) :: this
      real(dp), intent(in) :: omega

      associate( &
         nx => this%grid%nx, ny => this%grid%ny, &
         grid => this%grid)

#if defined(__GFORTRAN__) || defined(__NVCOMPILER)
            call lbm_collide_and_stream_fused(nx,ny,fsrc=grid%f1,fdst=grid%f2,omega=omega)
            call lbm_periodic_bc_push(nx,ny,grid%f2)
            swap: block
               real(dp), allocatable :: ftmp(:,:,:)
               call move_alloc(from=grid%f1,to=ftmp)
               call move_alloc(from=grid%f2,to=grid%f1)
               call move_alloc(from=ftmp,to=grid%f2)
            end block swap
#else
         select case(this%flip)
         case(0)
            call lbm_collide_and_stream_fused(nx,ny,fsrc=grid%f1,fdst=grid%f2,omega=omega)
            call lbm_periodic_bc_push(nx,ny,grid%f2)
         case(1)
            call lbm_collide_and_stream_fused(nx,ny,fsrc=grid%f2,fdst=grid%f1,omega=omega)
            call lbm_periodic_bc_push(nx,ny,grid%f1)
         end select
         this%flip = 1 - this%flip
#endif
      end associate

   end subroutine

   subroutine slbm_vars(this,rho,u)
      class(sim_slbm), intent(in) :: this
      real(dp), intent(out), contiguous :: rho(:,:), u(:,:,:)
      associate(grid => this%grid)
#if defined(__GFORTRAN__) || defined(__NVCOMPILER)
         call lbm_macros(grid%nx,grid%ny,1,grid%f1,rho,u(:,:,1),u(:,:,2))
#else
         select case(this%flip)
         case(0)
            call lbm_macros(grid%nx,grid%ny,1,grid%f1,rho,u(:,:,1),u(:,:,2))
         case(1)
            call lbm_macros(grid%nx,grid%ny,1,grid%f2,rho,u(:,:,1),u(:,:,2))
         end select
#endif
      end associate
   end subroutine

!
! PLUGIN INTERFACE
!
!void *siminit(
!   int nx, int ny, double dt,
!   double *rho, double *u, double *sigma, void *params);
function c_slbm_init(nx,ny,dt,rho,u,sigma,params) bind(c)
   integer(c_int), value :: nx, ny
   real(c_double), value :: dt
   real(c_double), intent(in) :: rho(nx,ny), u(nx,ny,2)
   real(c_double), intent(in), optional :: sigma(nx,ny,3)
   type(c_ptr), value :: params
   type(c_ptr) :: c_slbm_init

   type(sim_slbm), pointer :: sim_ptr 
   
   allocate(sim_ptr)
   call sim_ptr%init(nx,ny,dt,rho,u,sigma)
   c_slbm_init = c_loc(sim_ptr)

end function

!void simstep(void *sim, double omega);
subroutine c_slbm_step(ptr,omega) bind(c)
   type(c_ptr), value :: ptr
   real(c_double), value :: omega

   type(sim_slbm), pointer :: sim_ptr
   call c_f_pointer(ptr, sim_ptr)
   
   call sim_ptr%step(omega)
end subroutine

!void simvars(void *sim, double *rho, double *u);
subroutine c_slbm_vars(ptr,rho,u) bind(c)
   type(c_ptr), value :: ptr
   real(c_double), intent(out), target :: rho(*), u(*)

   real(c_double), pointer, contiguous :: rho_(:,:), u_(:,:,:)
   type(sim_slbm), pointer :: sim_ptr 

   call c_f_pointer(ptr, sim_ptr)
   
   associate(nx => sim_ptr%grid%nx, ny => sim_ptr%grid%ny)

   rho_(1:nx,1:ny) => rho(1:nx*ny)
   u_(1:nx,1:ny,1:2) => u(1:nx*ny*2)
   call sim_ptr%vars(rho_,u_)

   end associate

end subroutine


!void simfree(void *sim);
subroutine c_slbm_free(ptr) bind(c)
   type(c_ptr), value :: ptr
   
   !TODO: verify if this can be polymorphic
   type(sim_slbm), pointer :: simp 
   
   if(c_associated(ptr)) then
      call c_f_pointer(ptr, simp)
      DEALLOCATE(simp)
   end if

end subroutine


   function c_slbm_norm(nx,ny,u,ua) bind(c)
      integer(c_int), value :: nx, ny
      real(c_double), intent(in) :: u(nx,ny), ua(nx,ny)
      real(c_double) :: c_slbm_norm

      c_slbm_norm = norm2(u - ua) / norm2(ua)

   end function

end module
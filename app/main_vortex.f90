program main

   use fvm_bardow

   use barotropic_vortex_case, only: vortex_case_t

   use collision_bgk, only: collide_bgk
   use collision_trt, only: collide_trt

   ! use collision_bgk_improved, only: collide_bgk_improved

   use collision_regularized, only: collide_rr
   use periodic_lbm, only: perform_lbm_step, lbm_stream

   !$ use omp_lib, only: omp_get_wtime

   implicit none

   integer, parameter :: nx = 128, ny = 128
   integer, parameter :: nprint = 10000
   integer :: step, nsteps

   type(vortex_case_t) :: case
   real(wp) :: xc, yc, Rc, kappa

   type(lattice_grid) :: grid
   real(wp) :: cfl, dt, nu, tau, U0
   real(wp) :: Re

   real(wp) :: t, tmax

   real(wp) :: dt_over_tau
   character(len=64) :: arg

   integer, parameter :: dp = kind(1.0d0)
   real(dp) :: sbegin, send

   call alloc_grid(grid, nx, ny, nf = 2)

   grid%filename = "results"
   call grid%set_output_folder(foldername="vortex")

   grid%collision => collide_bgk
   grid%streaming => stream_fvm_bardow

   grid%logger => my_logger

   ! Advective Mach number
   ! umax = Mach * cs
   U0 = 0.1_wp / sqrt(3._wp)

   ! Vortex Mach number
   kappa = 0.2_wp / sqrt(3._wp)

   nu = 0.00001_wp

   Re = kappa*real(nx,wp)/nu
   print *, "Reynolds = ", Re

   ! tau = nu / cs**2, cs**2 = 1/3
   tau = 3.0_wp * nu 
   tau = 0.0001_wp

   tmax = (real(nx, wp) / U0) * 4.0_wp

   ! read value for dt/tau or cfl
   call get_command_argument(1,arg)
   
   read(arg,*) cfl
   dt = cfl
   dt_over_tau = dt/tau

   !read(arg,*) dt_over_tau
   !dt = dt_over_tau*tau
   !cfl = dt
   
   print *, "tau = ", tau
   print *, "dt/tau = ", dt/tau
   print *, "cfl = ", cfl

   call set_properties(grid, nu, dt, magic=1._wp/4._wp)

   print *, "omega = ", grid%omega

   ! ---- prepare flow case ----

   xc = real(nx,wp)/2._wp
   yc = real(ny,wp)/2._wp

   Rc = real(nx,wp)/10._wp

   case = vortex_case_t(U0,xc,yc,Rc, &
      eps = kappa)
   
   ! ---- number of time steps

   nsteps = int(tmax / dt) + 1
   print*, "nsteps = ", nsteps

   t = 0._wp

   call apply_initial_condition(case, grid)

   call output_gnuplot(grid,step=0)
   call output_npy(grid,step=0)
   call output_vtk(grid,step=0)

   call grid%logger(step=0)

   !$ sbegin = omp_get_wtime()

   time: do step = 1, nsteps

      call perform_step(grid)
      !call perform_lbm_step(grid)
      t = t + dt

      if (mod(step,nprint) == 0) then
         print *, "step = ", step
         call update_macros(grid)
         call output_gnuplot(grid,step)
         call output_npy(grid,step)
         call output_vtk(grid,step)

         call grid%logger(step)
      end if

      if (t >= tmax) then
         exit time
      end if

   end do time

   !$ send = omp_get_wtime()
   !$ print *, "MLUPS ", (real(nx,wp) * real(ny,wp) * real(step,wp) * 1.e-6_wp) / (send - sbegin)

   ! calculate average L2-norm
   print *, "Final time = ", t

   call dealloc_grid(grid)

contains

   subroutine apply_initial_condition(case, grid)
      type(vortex_case_t), intent(in) :: case
      type(lattice_grid), intent(inout) :: grid

      real(wp), parameter :: rho0 = 1.0_wp

      call case%eval(grid%nx,grid%ny,grid%rho,grid%ux,grid%uy)

      call set_pdf_to_equilibrium(grid)

   end subroutine

   subroutine my_logger(grid,step)
      class(lattice_grid), intent(in) :: grid
      integer, intent(in) :: step

      write(grid%logunit, *) step, step*dt, maxval(hypot(grid%ux,grid%uy))
      flush(grid%logunit)

   end subroutine

end program
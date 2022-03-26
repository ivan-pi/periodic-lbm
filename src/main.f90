program main

   use fvm_bardow
   use taylor_green, only: taylor_green_t, pi
   use collision_trt, only: collide_trt
   use collision_bgk_improved, only: collide_bgk_improved

   !$ use omp_lib, only: omp_get_wtime

   implicit none

   integer, parameter :: nx = 256, ny = 256
   integer, parameter :: nprint = 20000
   integer :: step, nsteps

   type(taylor_green_t) :: tg

   type(lattice_grid) :: grid
   real(wp) :: cfl, dt, nu, tau
   real(wp) :: kx, ky, umax, nrm

   real(wp) :: t, tmax

   real(wp) :: dt_over_tau
   character(len=64) :: arg

   integer, parameter :: dp = kind(1.0d0)
   real(dp) :: sbegin, send

   call alloc_grid(grid, nx, ny)

   grid%filename = "results"
   call grid%set_output_folder(foldername="taylor_green")

   grid%collision => collide_bgk
   grid%streaming => stream_fdm_bardow

   grid%logger => my_logger

   ! umax = Mach * cs
   umax = 0.01_wp / sqrt(3._wp)

   ! nu = (umax * L) / Re
   nu = (umax * real(nx,wp)) / 100._wp
   
   ! tau = nu / cs**2, cs**2 = 1/3
   tau = 3.0_wp * nu 

   ! read value for dt/tau
   call get_command_argument(1,arg)
   read(arg,*) cfl
   dt = cfl/sqrt(2.0_wp)
   dt_over_tau = dt/tau

   !read(arg,*) dt_over_tau
   !dt = dt_over_tau*tau
   !cfl = 0.1_wp
   !cfl = sqrt(2._wp)*dt
   
   print *, "tau = ", tau
   print *, "dt/tau = ", dt/tau
   print *, "cfl = ", cfl

   call set_properties(grid, nu, dt, magic=1._wp/4._wp)

   print *, "omega = ", grid%omega

   ! ---- prepare flow case ----

   kx = 2*pi/real(nx,wp)
   ky = 2*pi/real(ny,wp)

   tg = taylor_green_t(nx,ny,kx,ky,umax,nu)
   print *, "umax = ", umax
   print *, "tc   = ", tg%td
   call write_gnuplot_include()

   tmax = log(2._wp)*tg%decay_time()
   nsteps = int(1.1_wp*tmax/dt)
   print*, "nsteps = ", nsteps

   t = 0._wp
   call apply_initial_condition(tg, grid)

   call output_grid_txt(grid,step=0)
   call output_vtk(grid,step=0)
   
   call grid%logger(step=0)

   !$ sbegin = omp_get_wtime()

   time: do step = 1, nsteps

      call perform_step(grid)
      t = t + dt

      if (mod(step,nprint) == 0) then
         print *, step
         call update_macros(grid)
         call output_grid_txt(grid,step)
         call output_vtk(grid,step)
         call grid%logger(step)
      end if

      if (t >= tmax) then
         call update_macros(grid)
         call output_grid_txt(grid,step)
         call output_vtk(grid,step)
         call grid%logger(step)
         exit time
      end if
   end do time

   !$ send = omp_get_wtime()
   !$ print *, "MLUPS ", (real(nx,wp) * real(ny,wp) * real(step,wp) * 1.e-6_wp) / (send - sbegin)

   ! calculate average L2-norm
   nrm = calc_L2_norm(tg, grid, t)
   print *, "L2-norm = ", nrm
   print *, "Final time = ", t

   call dealloc_grid(grid)

contains

   subroutine apply_initial_condition(case, grid)
      type(taylor_green_t), intent(in) :: case
      type(lattice_grid), intent(inout) :: grid

      real(wp), parameter :: rho0 = 1.0_wp

      call case%eval(t=0.0_wp, &
                   p=grid%rho, &
                   ux=grid%ux, &
                   uy=grid%uy)

      ! convert pressure to lattice density
      grid%rho = grid%rho/grid%csqr + rho0

      call set_pdf_to_equilibrium(grid)

   end subroutine

   subroutine my_logger(grid,step)
      class(lattice_grid), intent(in) :: grid
      integer, intent(in) :: step

      write(grid%logunit, *) step, step*dt, maxval(hypot(grid%ux,grid%uy))
      flush(grid%logunit)

   end subroutine

   subroutine write_gnuplot_include()

      integer :: unit

      open(newunit=unit,file="lattice_grid_log.incl",status='unknown')

      write(unit,*) "dt = ", dt
      write(unit,*) "umax = ", umax
      write(unit,*) "tc = ", tg%td

      close(unit)
   end subroutine


   function calc_L2_norm(case, grid, t) result(nrm)
      type(taylor_green_t), intent(in) :: case
      type(lattice_grid), intent(in) :: grid
      real(wp), intent(in) :: t
      real(wp) :: nrm

      real(wp), allocatable :: pa(:,:), uxa(:,:), uya(:,:)
      real(wp) :: above, below

      allocate(pa, mold=grid%rho) ! not needed
      allocate(uxa, mold=grid%ux)
      allocate(uya, mold=grid%uy)

      call case%eval(t=t, &
                   p=pa, &
                   ux=uxa, &
                   uy=uya)

      above = norm2(hypot(grid%ux-uxa, grid%uy-uya))
      below = norm2(hypot(uxa, uya))
      nrm = above/below

      ! # Code used in the Python version
      ! def L2_error(u,v,ua,va):
      !     return np.sqrt(np.sum((u-ua)**2 + (v-va)**2)/np.sum(ua**2 + va**2))

      !above = sum((grid%ux - uxa)**2 + (grid%uy - uya)**2)
      !below = sum(uxa**2 + uya**2)
      !nrm = sqrt(above/below)

   end function

end program
module collision_regularized

   use precision, only: wp
   use fvm_bardow, only: lattice_grid

   implicit none
   private

   public :: collide_rr

   real(wp), parameter :: w0 = 4.0_wp/9.0_wp
   real(wp), parameter :: ws = 1.0_wp/9.0_wp
   real(wp), parameter :: wd = 1.0_wp/36.0_wp
   real(wp), parameter :: csqr = 1.0_wp/3.0_wp

contains

   ! Recursively regularized collision
   !
   !
   subroutine collide_rr(grid)
      class(lattice_grid), intent(inout) :: grid

      integer :: ld

      ld = size(grid%f,1)

!#if SPLIT
!      call rr_kernel_split(grid%nx, grid%ny, ld, &
!         grid%f(:,:,:,grid%inew), &
!         grid%omega)
!#else
      call rr_kernel_naive(grid%nx, grid%ny, ld, &
         grid%f(:,:,:,grid%inew), &
         grid%omega)
!#endif

   contains

      subroutine rr_kernel_naive(nx,ny,ld,f1,omega)
         integer, intent(in) :: nx, ny, ld
         real(wp), intent(inout) :: f1(ld,nx,0:8)
         real(wp), intent(in) :: omega

         real(wp) :: rho, invrho, ux, uy

         real(wp) :: omega_w0, omega_ws, omega_wd
         real(wp) :: vC, vE, vN, vW, vS, vNE, vNW, vSW, vSE
         
         real(wp) :: feq(0:8)
         real(wp) :: uxx, uyy, uxxy, uyyx, uxxyy
         real(wp) :: uxpy, uxmy, u3p, u3m

         real(wp) :: indp0, indps, indpd, indp57, indp68

         real(wp) :: axx, axy, ayy, axxy, ayyx, axxyy, tmp, a3p, a3m


         integer :: x, y

      !$omp parallel default(private) shared(f1,omega,nx,ny)

         omega_w0 = w0*(1.0_wp - omega) 
         omega_ws = ws*(1.0_wp - omega) 
         omega_wd = wd*(1.0_wp - omega) 

         !$omp do schedule(static)
         do x = 1, nx 
            do y = 1, ny

            ! pull pdfs
            vC  = f1(y,x,0)
            VE  = f1(y,x,1)
            vN  = f1(y,x,2)
            vW  = f1(y,x,3)
            vS  = f1(y,x,4)
            vNE = f1(y,x,5)
            vNW = f1(y,x,6)
            vSW = f1(y,x,7)
            vSE = f1(y,x,8)

            !
            ! macroscopic values
            !
            rho = (((vNE + vSW) + (vNW + vSE)) + &
                   ((vE + vW) + (vN + vS))) + vC

            invrho = 1.0_wp / rho

            ux = invrho * (((vNE - vSW) + (vSE - vNW)) + (vE - vW))
            uy = invrho * (((vNE - vSW) + (vNW - vSE)) + (vN - vS))

            uxx   = ux*ux
            uyy   = uy*uy
            uxxy  = uxx*uy
            uyyx  = uyy*ux
            uxxyy = uxx*uyy

            indp0 = 1.0_wp - 1.5_wp * (uxx + uyy)
            indps = indp0 - 4.5_wp*uxxyy
            indpd = indp0 + 9.0_wp*uxxyy

            indp0 = indp0 + 2.25_wp*uxxyy

            !
            ! equilibrium parts, 0 - 4
            !
            feq(0) = w0*rho*indp0

            feq(1) = ws*rho*(indps + 3.0_wp*ux + 4.5_wp*(uxx - uyyx))
            feq(3) = ws*rho*(indps - 3.0_wp*ux + 4.5_wp*(uxx + uyyx))

            feq(2) = ws*rho*(indps + 3.0_wp*uy + 4.5_wp*(uyy - uxxy))
            feq(4) = ws*rho*(indps - 3.0_wp*uy + 4.5_wp*(uyy + uxxy))
            
            !
            ! non-equilibrium parts, 0 - 4
            !
            vC = vC - feq(0)
            vE = vE - feq(1)
            vN = vN - feq(2)
            vW = vW - feq(3)
            vS = vS - feq(4)

            axx = csqr*(2*(vE + vW) - (vN + vS) - vC)
            ayy = csqr*(2*(vN + vS) - (vE + vW) - vC)

            !
            ! equilibrium parts, 5 - 8
            !
            u3p = uxxy + uyyx
            uxpy = ux + uy
            indp57 = indpd + 4.5_wp*uxpy*uxpy
            
            feq(5) = wd*rho*(indp57 + 3.0_wp*uxpy + 9.0_wp*u3p)
            feq(7) = wd*rho*(indp57 - 3.0_wp*uxpy - 9.0_wp*u3p)

            u3m = uxxy - uyyx
            uxmy = ux - uy
            indp68 = indpd + 4.5_wp*uxmy*uxmy

            feq(6) = wd*rho*(indp68 - 3.0_wp*uxmy + 9.0_wp*u3m)
            feq(8) = wd*rho*(indp68 + 3.0_wp*uxmy - 9.0_wp*u3m)

            ! ----------------------

            !
            ! non-equilibrium parts, 5 - 8
            !
            vNE = vNE - feq(5) 
            vNW = vNW - feq(6)
            vSW = vSW - feq(7)
            vSE = vSE - feq(8)

            tmp = 2.0_wp*csqr*(vNE + vNW + vSW + vSE)
            
            axx = axx + tmp
            ayy = ayy + tmp

            axy = ((vNE + vSW) - (vNW + vSE))

            axxy = 2.0_wp*ux*axy + uy*axx
            ayyx = 2.0_wp*uy*axy + ux*ayy
            axxyy = 2.0_wp*(ux*ayyx + uy*axxy) - uxx*ayy - uyy*axx -  4.0_wp*ux*uy*axy

            indp0 = -1.5_wp*(axx + ayy)
            indps = indp0 - 4.5_wp*axxyy 
            indpd = 9.0_wp*axxyy - 2.0_wp*indp0

            indp0 = indp0 + 2.25_wp*axxyy

            vC = indp0

            vE = indps + 4.5_wp*(axx - ayyx)
            vW = indps + 4.5_wp*(axx + ayyx)
            
            vN = indps + 4.5_wp*(ayy - axxy)
            vS = indps + 4.5_wp*(ayy + axxy)
            
            vNE = indpd + 9.0_wp*(axxy + ayyx + axy)
            vSW = indpd - 9.0_wp*(axxy + ayyx - axy)

            vNW = indpd + 9.0_wp*(axxy - ayyx - axy)
            vSE = indpd - 9.0_wp*(axxy - ayyx + axy)

            f1(y,x,0) = feq(0) + omega_w0 * vC
            f1(y,x,1) = feq(1) + omega_ws * vE
            f1(y,x,2) = feq(2) + omega_ws * vN
            f1(y,x,3) = feq(3) + omega_ws * vW
            f1(y,x,4) = feq(4) + omega_ws * vS
            f1(y,x,5) = feq(5) + omega_wd * vNE
            f1(y,x,6) = feq(6) + omega_wd * vNW
            f1(y,x,7) = feq(7) + omega_wd * vSW
            f1(y,x,8) = feq(8) + omega_wd * vSE

            end do
         end do
         !$omp end do

      !$omp end parallel

      end subroutine


      ! Turns out the split version is slower than the naive one
      ! probably due to reading each pdf 3 times, instead of only 2 
      ! in the naive version
      !
      ! TODO: some mistakes remain in this version. Use the naive version.
      !
      subroutine rr_kernel_split(nx,ny,ld,f1,omega)
         integer, intent(in) :: nx, ny, ld
         real(wp), intent(inout) :: f1(ld,nx,0:8)
         real(wp), intent(in) :: omega

         real(wp) :: rho(ny), invrho, ux(ny), uy(ny)

         real(wp) :: omega_w0, omega_ws, omega_wd
         real(wp) :: vC, vE, vN, vW, vS, vNE, vNW, vSW, vSE
         
         real(wp) :: uxx, uyy, uxxy, uyyx, uxxyy(ny)
         real(wp) :: uxpy, uxmy, u3p, u3m

         real(wp) :: indp0(ny)
         real(wp) :: indps, indpd

         real(wp) :: axx(ny), axy(ny), ayy(ny)
         real(wp) :: axxy(ny), ayyx(ny), axxyy(ny)
         real(wp) :: tmp, a3p, a3m

         integer :: x, y

      !$omp parallel default(private) shared(f1,omega,nx,ny)

         omega_w0 = w0*(1.0_wp - omega) 
         omega_ws = ws*(1.0_wp - omega) 
         omega_wd = wd*(1.0_wp - omega) 

         !$omp do schedule(static)
         do x = 1, nx  
            do y = 1, ny

               ! pull pdfs
               vC  = f1(y,x,0)
               VE  = f1(y,x,1)
               vN  = f1(y,x,2)
               vW  = f1(y,x,3)
               vS  = f1(y,x,4)
               vNE = f1(y,x,5)
               vNW = f1(y,x,6)
               vSW = f1(y,x,7)
               vSE = f1(y,x,8)

               !
               ! macroscopic values
               !
               rho(y) = (((vNE + vSW) + (vNW + vSE)) + &
                      ((vE + vW) + (vN + vS))) + vC

               invrho = 1.0_wp / rho(y)

               ux(y) = invrho * (((vNE - vSW) + (vSE - vNW)) + (vE - vW))
               uy(y) = invrho * (((vNE - vSW) + (vNW - vSE)) + (vN - vS))

               uxx = ux(y)*ux(y)
               uyy = uy(y)*uy(y)
               uxxyy(y) = uxx*uyy

               indp0(y) = 1.0_wp - 1.5_wp * (uxx + uyy)
               f1(y,x,0) = w0*rho(y)*(indp0(y) - 2.25_wp*uxxyy(y))

               axx(y) = f1(y,x,0) - vC 
               ayy(y) = axx(y)
            end do


            do y = 1, ny

               indps = indp0(y) - 4.5_wp*uxxyy(y)
               uxx = ux(y)*ux(y)
               uyyx = uy(y)*uy(y)*ux(y)

               vE = f1(y,x,1) 
               vW = f1(y,x,3)
               f1(y,x,1) = ws*rho(y)*(indps + 3.0_wp*ux(y) + 4.5_wp*(uxx - uyyx))
               f1(y,x,3) = ws*rho(y)*(indps - 3.0_wp*ux(y) + 4.5_wp*(uxx + uyyx))

               tmp = (vE - f1(y,x,1)) + (vW - f1(y,x,3))
               axx(y) = axx(y) + 2*tmp
               ayy(y) = ayy(y) - tmp

            end do

            do y = 1, ny

               indps = indp0(y) - 4.5_wp*uxxyy(y)
               uyy = uy(y)*uy(y)
               uxxy = ux(y)*ux(y)*uy(y)

               vN = f1(y,x,2) 
               vS = f1(y,x,4)
               f1(y,x,2) = ws*rho(y)*(indps + 3.0_wp*uy(y) + 4.5_wp*(uyy - uxxy))
               f1(y,x,4) = ws*rho(y)*(indps - 3.0_wp*uy(y) + 4.5_wp*(uyy + uxxy))

               tmp = (vN - f1(y,x,2)) + (vS - f1(y,x,4))
               axx(y) = axx(y) - tmp
               ayy(y) = ayy(y) + 2*tmp

            end do
            

            !
            ! equilibrium parts, 5 - 8
            !
            do y = 1, ny

               uxxy = ux(y)*ux(y)*uy(y)
               uyyx = uy(y)*uy(y)*ux(y)

               u3p = uxxy + uyyx
               uxpy = ux(y) + uy(y)
               indpd = indp0(y) + 9.0_wp*uxxyy(y) + 4.5_wp*uxpy*uxpy
            
               vNE = f1(y,x,5)
               vSW = f1(y,x,7)
               f1(y,x,5) = wd*rho(y)*(indpd + 3.0_wp*uxpy + 9.0_wp*u3p)
               f1(y,x,7) = wd*rho(y)*(indpd - 3.0_wp*uxpy - 9.0_wp*u3p)

               vNE = vNE - f1(y,x,5)
               vSW = vSW - f1(y,x,7)

               tmp = vNE + vSW
               axx(y) = axx(y) + 2*tmp
               axy(y) = tmp
               ayy(y) = ayy(y) + 2*tmp
            end do

            do y = 1, ny

               uxxy = ux(y)*ux(y)*uy(y)
               uyyx = uy(y)*uy(y)*ux(y)
               
               u3m = uxxy - uyyx
               uxmy = ux(y) - uy(y)
               indpd = indp0(y) + 9.0_wp*uxxyy(y) + 4.5_wp*uxmy*uxmy

               vNW = f1(y,x,6)
               vSE = f1(y,x,8)
               f1(y,x,6) = wd*rho(y)*(indpd - 3.0_wp*uxmy + 9.0_wp*u3m)
               f1(y,x,8) = wd*rho(y)*(indpd + 3.0_wp*uxmy - 9.0_wp*u3m)

               vNW = vNW - f1(y,x,6)
               vSE = vSE - f1(y,x,8)

               tmp = vNW + vSE
               axx(y) = axx(y) + 2*tmp
               axy(y) = axy(y) - tmp
               ayy(y) = ayy(y) + 2*tmp

               axx(y) = csqr*axx(y)
               ayy(y) = csqr*ayy(y)
            end do

            !
            ! non-equilibrium parts
            !

            do y = 1, ny

               uxx = ux(y)*ux(y)
               uyy = uy(y)*uy(y)

               axxy(y) = 2.0_wp*ux(y)*axy(y) + uy(y)*axx(y)
               ayyx(y) = 2.0_wp*uy(y)*axy(y) + ux(y)*ayy(y)

               axxyy(y) = 2.0_wp*(ux(y)*ayyx(y) + uy(y)*axxy(y)) - &
                  uxx*ayy(y) - uyy*axx(y) -  4.0_wp*ux(y)*uy(y)*axy(y)

               indp0(y) = -1.5_wp*(axx(y) + ayy(y))

               vC = indp0(y) + 2.25_wp*axxyy(y)
               f1(y,x,0) = f1(y,x,0) + omega_w0*vC

            end do

            do y = 1, ny
               indps = indp0(y) - 4.5_wp*axxyy(y) 
               vE = indps + 4.5_wp*(axx(y) - ayyx(y))
               vW = indps + 4.5_wp*(axx(y) + ayyx(y))
               f1(y,x,1) = f1(y,x,1) + omega_ws * vE
               f1(y,x,3) = f1(y,x,3) + omega_ws * vW
            end do

            do y = 1, ny
               indps = indp0(y) - 4.5_wp*axxyy(y)
               vN = indps + 4.5_wp*(ayy(y) - axxy(y))
               vS = indps + 4.5_wp*(ayy(y) + axxy(y))
               f1(y,x,2) = f1(y,x,2) + omega_ws * vN
               f1(y,x,4) = f1(y,x,4) + omega_ws * vS
            end do

            do y = 1, ny
               indpd = 9.0_wp*axxyy(y) - 2.0_wp*indp0(y)
               a3p = axxy(y) + ayyx(y)
               vNE = indpd + 9.0_wp*(a3p + axy(y))
               vSW = indpd - 9.0_wp*(a3p + axy(y))
               f1(y,x,5) = f1(y,x,5) + omega_wd * vNE 
               f1(y,x,7) = f1(y,x,7) + omega_wd * vSW
            end do

            do y = 1, ny
               indpd = 9.0_wp*axxyy(y) - 2.0_wp*indp0(y)
               a3m = axxy(y) - ayyx(y)
               vNW = indpd + 9.0_wp*(a3m + axy(y))
               vSE = indpd - 9.0_wp*(a3m + axy(y))
               f1(y,x,6) = f1(y,x,6) + omega_wd * vNW
               f1(y,x,8) = f1(y,x,8) + omega_wd * vSE
            end do

         end do
         !$omp end do

      !$omp end parallel

      end subroutine

   end subroutine






end module

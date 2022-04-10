module vorticity

    use precision, only: wp

    implicit none
    private

    public :: vorticity_2nd
    public :: vorticity_4th

contains

    subroutine vorticity_2nd(ux,uy,omega)
        real(wp), intent(in) :: ux(:,:), uy(:,:)
        real(wp), intent(out) :: omega(:,:)

        real(wp) :: duydx, duxdy
        integer :: xp1, xm1, yp1, ym1
        integer :: x, y, nx, ny

        ny = size(ux,1)
        nx = size(ux,2)

        !$omp parallel default(private) shared(ux,uy,omega,nx,ny)
        !$omp do schedule(static) 
        do x = 1, nx
            xp1 = mod(x,nx) + 1
            xm1 = mod(nx + x - 2,nx) + 1
            do y = 1, ny
                yp1 = mod(y,ny) + 1
                ym1 = mod(ny + y - 2,ny) + 1

                duydx = 0.5_wp*(uy(y,xp1) - uy(y,xm1))
                duxdy = 0.5_wp*(ux(yp1,x) - ux(ym1,x))

                omega(y,x) = duydx - duxdy

            end do
        end do
        !$omp end do 
        !$omp end parallel

    end subroutine


    subroutine vorticity_4th(ux,uy,omega)
        real(wp), intent(in) :: ux(:,:), uy(:,:)
        real(wp), intent(out) :: omega(:,:)

        real(wp) :: duydx, duxdy
        integer :: xp1, xm1, yp1, ym1, xp2, yp2, xm2, ym2
        integer :: x, y

        real(wp), parameter :: t1 = 1._wp/12._wp, &
                               t2 = 2._wp/3._wp
        integer :: nx, ny

        nx = size(ux,2)
        ny = size(ux,1)

        !$omp parallel default(private) shared(ux,uy,omega,nx,ny)
        !$omp do schedule(static) 
        do x = 1, nx
            
            xp1 = mod(x,nx) + 1
            xm1 = mod(nx + x - 2,nx) + 1
            xp2 = mod(x + 1,nx) + 1
            xm2 = mod(nx + x - 3,nx) + 1

            do y = 1, ny

                yp1 = mod(y,ny) + 1
                ym1 = mod(ny + y - 2,ny) + 1
                yp2 = mod(y + 1,ny) + 1
                ym2 = mod(ny + y - 3,ny) + 1

                duydx = t1*(uy(y,xp1) - uy(y,xm1)) + t2*(uy(y,xm2) - uy(y,xp2))
                duxdy = t1*(ux(yp1,x) - ux(ym1,x)) + t2*(ux(ym2,x) - ux(yp2,x))

                omega(y,x) = duydx - duxdy

            end do
        end do
        !$omp end do 
        !$omp end parallel

    end subroutine
    
end module
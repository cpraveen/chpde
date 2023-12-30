subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    ! Set initial conditions for q.

    implicit none
    include 'param.h'

    integer, intent(in) :: meqn, mbc, mx, my, maux
    real(kind=8), intent(in) :: xlower, ylower, dx, dy
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

    integer :: i, j

    do j=1,my
        do i=1,mx
            q(1,i,j) = rhoinf
            q(2,i,j) = rhoinf * uinf
            q(3,i,j) = rhoinf * vinf
            q(4,i,j) = pinf/(gamma - 1.0d0) &
                       + 0.5d0 * rhoinf * (uinf**2 + vinf**2)
        end do
    end do
end subroutine qinit

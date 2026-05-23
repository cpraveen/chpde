subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.

    implicit none
    
    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    integer :: i
    real(kind=8) :: xcell
    real(kind=8) :: rhol, vl, pl, gammal
    real(kind=8) :: rhor, vr, pr, gammar
    common /cparam/ rhol, vl, pl, gammal, rhor, vr, pr, gammar

    do i=1,mx
        xcell = xlower + (i-0.5d0)*dx
        if (xcell .lt. 0.5d0) then
           q(1,i) = rhol
           q(2,i) = rhol * vl
           q(3,i) = pl/(gammal - 1.0d0) + 0.5d0 * rhol * vl**2
           q(4,i) = 1.0d0/(gammal - 1.0d0)
        else
           q(1,i) = rhor
           q(2,i) = rhor * vr
           q(3,i) = pr/(gammar - 1.0d0) + 0.5d0 * rhor * vr**2
           q(4,i) = 1.0d0/(gammar - 1.0d0)
        endif
    enddo

end subroutine qinit


subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

   ! Set initial conditions for q.

   implicit none
   include 'param.h'

   integer, intent(in) :: meqn, mbc, mx, my, maux
   real(kind=8), intent(in) :: xlower, ylower, dx, dy
   real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
   real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

   integer :: i, j, m

   do j=1,my
      do i=1,mx
         do m=1,meqn
            q(m,i,j) = qinf(m)
         enddo
      enddo
   enddo

end subroutine qinit

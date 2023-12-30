subroutine mapc2p(xc, yc, xp, yp)
   implicit none
   include 'param.h'
   real(8),intent(in)    :: xc, yc
   real(8),intent(inout) :: xp, yp
   ! Local variables
   real(8) :: r, theta

   r     = rin + xc * (rout - rin)
   theta = 2.0d0 * M_PI * yc

   xp = r * dcos(theta)
   yp = r * dsin(theta)

end subroutine mapc2p

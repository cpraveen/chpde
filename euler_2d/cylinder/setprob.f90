subroutine setprob

   implicit none
   include 'param.h'

   M_PI = 4.0d0 * datan(1.0d0)

   gamma = 1.4d0

   machinf = 0.1d0

   rhoinf = 1.0d0
   uinf = 1.0d0
   vinf = 0.0d0
   pinf = 1.0d0 / (gamma * machinf**2)

   qinf(1) = rhoinf
   qinf(2) = rhoinf * uinf
   qinf(3) = rhoinf * vinf
   qinf(4) = pinf/(gamma-1.0d0) + 0.5d0 * rhoinf * (uinf**2 + vinf**2)

end subroutine setprob

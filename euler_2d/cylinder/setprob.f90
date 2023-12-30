subroutine setprob

   implicit none
   include 'param.h'

   M_PI = 4.0d0 * atan(1.0d0)

   gamma = 1.4d0
   machinf = 0.1d0

   rin = 1.0d0
   rout = 50.0d0 * rin

   rhoinf = 1.0d0
   uinf = 1.0d0
   vinf = 0.0d0
   pinf = 1.0d0 / (gamma * machinf**2)

end subroutine setprob

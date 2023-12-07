      subroutine setprob
      implicit double precision (a-h,o-z)
      common /comtraf/ u1,u2,rho1,rho2
c
c     # Set the velocity for scalar advection
c     # This value is passed to the Riemann solver rp1.f in a common block
c
c
      call opendatafile(7,'setprob.data')

      read(7,*) u1,u2
      read(7,*) rho1,rho2

      return
      end

      program driver
c
c  Generic driver routine for claw1ez
c
c  Author: Randall J. LeVeque
c  Version of March, 1999 --  CLAWPACK Version 4.0
c
c
      implicit double precision (a-h,o-z)

c     # set parameters for maximum array sizes used in declarations
c     # these must be increased for larger problems.
c
c
      parameter (maxmx = 4000)
      parameter (mwork = 152104)

      parameter (mbc = 2)
      parameter (meqn = 3)
      parameter (mwaves = 3)
      parameter (maux = 0)
c       # NOTE: if maux>0 you must declare aux properly below!
c
      dimension q(1-mbc:maxmx+mbc, meqn)

      dimension  aux(1)   !# dummy variable since no aux arrays used
c     dimension  aux(1-mbc:maxmx+mbc, maux)

      dimension work(mwork)
      dimension mthlim(mwaves)

      call claw1ez(maxmx,meqn,mwaves,mbc,maux,mwork,mthlim,
     &           q,work,aux)

      stop 
      end

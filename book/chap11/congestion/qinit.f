c
c
c =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)
      common /comic/ beta
c
c
      do 150 i=1,mx
	 xcell = xlower + (i-0.5d0)*dx
	 q(i,1) = 0.25d0 + 0.7d0*dexp(-beta * (xcell-0.0d0)**2)
  150    continue
c
      return
      end

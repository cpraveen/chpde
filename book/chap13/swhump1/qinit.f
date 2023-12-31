c
c
c =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c     # Pulse in pressure, zero velocity
c
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)
      common /cqinit/ beta,ic
c
c
      do 150 i=1,mx
	 xcell = xlower + (i-0.5d0)*dx

	 go to (10,20,30) ic


   10    continue
c        # half ellipse:
	 if (xcell.gt.-4d0 .and. xcell.lt.-2d0) then
	     q(i,1) = 1.d0 + dsqrt(1.d0 - (xcell+3.d0)**2)
	    else
	     q(i,1) = 1.d0
	    endif
	 q(i,2) = 0.d0
         go to 150

   20    continue
c        # single discontinuity:
         if (xcell .lt. -2.d0) then
	     q(i,1) = 2.d0
	    else
	     q(i,1) = 1.d0
	    endif
	 q(i,2) = 0.d0
         go to 150

   30    continue
c        # Gaussian:
	 q(i,1) = 1.d0 + 0.4d0 * dexp(-beta*(xcell)**2)  
	 q(i,2) = 0.d0
         go to 150

  150    continue
c
      return
      end

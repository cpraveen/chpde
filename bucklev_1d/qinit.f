c
c
c =========================================================
       subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c
c
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
      
      common /comic/ ur, ul
c
c
      x0 = 0.0d0
c
c     # set Riemann initial values
      do i=1,mx
         xcell = xlower + (i-0.5d0)*dx
         if(xcell.lt.x0) then
            q(1,i) = ul
         else
            q(1,i)= ur
         endif
      enddo

c
      return
      end

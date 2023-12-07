c =========================================================
       subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
      common /comtraf/ u1,u2,rho1,rho2
c
      do i=1,mx
         xcell = xlower + (i-0.5d0)*dx
         if (xcell.lt.0.0d0) then
            q(1,i) = rho1
         else
            q(1,i) = rho2
         endif
      enddo
c
      return
      end

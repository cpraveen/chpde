c     ============================================
      subroutine setaux(mbc,mx,xlower,dx,maux,aux)
c     ============================================
c
c     # set auxiliary arrays 
c     # aux(1,i) = velocity u_i in i'th cell for advection
c
c     # periodic smooth variation u = sin(2*pi*x), evaluated
c     # at cell centers
c     
      implicit double precision (a-h,o-z)
      dimension aux(maux,1-mbc:mx+mbc)
      common /comtraf/ u1,u2,rho1,rho2

      do i=1-mbc,mx+mbc
         xedge = xlower + (i-1)*dx
         if (xedge .lt. 0.d0) then
            aux(i,1) = u1
         else
            aux(i,1) = u2
         endif
      enddo
c
      return
      end

c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                  dx,dy,q)
c     =====================================================
c
c     # Set initial conditions for q.
 
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
 
      do i=1,mx
         xi = xlower + (i-0.5d0)*dx
         do j=1,my
            yj = ylower + (j-0.5d0)*dy
            if (xi.lt.0.6d0 .and. xi.gt.0.1d0 .and. yj.gt.-0.25d0 .and.
     &          yj.lt.0.25d0) then
               q(1,i,j) = 1.d0
            else
               q(1,i,j) = 0.d0
            endif
            r = dsqrt((xi+0.45d0)**2 + yj**2)
            if (r .lt. 0.35d0) then
               q(1,i,j) = 1.d0 - r/0.35d0
            endif
         enddo
      enddo

      return
      end


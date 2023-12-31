c     ============================================
      subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays 
c     #   aux(1,i,j) is edge velocity at "left" boundary of grid point (i,j)
c     #   aux(2,i,j) is edge velocity at "bottom" boundary of grid point (i,j)
 
      implicit double precision (a-h,o-z)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      do i=1-mbc,mx+mbc
         do j=1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xll = xlower + (i-1)*dx
            yll = ylower + (j-1)*dy
c           # difference stream function psi to get normal velocities:
            aux(1,i,j) =  (psi(xll, yll+dy) - psi(xll,yll)) / dy
            aux(2,i,j) = -(psi(xll+dx, yll) - psi(xll,yll)) / dx
         enddo
      enddo
 
      return
      end


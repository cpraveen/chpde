subroutine setaux(mbc,mx,my,xlower,ylower,dxc,dyc,maux,aux)

   ! Called at start of computation before calling qinit, and
   ! when AMR is used, also called every time a new grid patch is created.
   ! Use to set auxiliary arrays aux(1:maux, 1-mbc:mx+mbc, 1-mbc:my+mbc).
   ! Note that ghost cell values may need to be set if the aux arrays
   ! are used by the Riemann solver(s).
   !
   !    aux(1,i,j)  = ax
   !    aux(2,i,j)  = ay   where (ax,ay) is unit normal to left face
   !    aux(3,i,j)  = ratio of length of left face to dyc

   !    aux(4,i,j)  = bx
   !    aux(5,i,j)  = by   where (bx,by) is unit normal to bottom face
   !    aux(6,i,j)  = ratio of length of bottom face to dxc

   !    aux(7,i,j)  = ratio of cell area to dxc*dyc
   !                  (approximately Jacobian of mapping function)
 
   implicit none
   integer, intent(in) :: mbc,mx,my,maux
   real(kind=8), intent(in) :: xlower,ylower,dxc,dyc
   real(kind=8), intent(out) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

   integer      :: i, j, ic
   real(kind=8) :: xccorn(5),yccorn(5),xpcorn(5),ypcorn(5)
   real(kind=8) :: ax, ay, bx, by, anorm, bnorm, area

   print*,'Setting bottom profile in aux ...'

   do j=1,my
      do i=1,mx
         ! computational points (xc,yc) are mapped to physical
         ! coordinates (xp,yp) by mapc2p:
 
         ! lower left corner:
         xccorn(1) = xlower + (i-1)*dxc
         yccorn(1) = ylower + (j-1)*dyc
         call mapc2p(xccorn(1),yccorn(1),xpcorn(1),ypcorn(1))

         ! upper left corner:
         xccorn(2) = xccorn(1)
         yccorn(2) = yccorn(1) + dyc
         call mapc2p(xccorn(2),yccorn(2),xpcorn(2),ypcorn(2))

         ! upper right corner:
         xccorn(3) = xccorn(1) + dxc
         yccorn(3) = yccorn(1) + dyc
         call mapc2p(xccorn(3),yccorn(3),xpcorn(3),ypcorn(3))

         ! lower right corner:
         xccorn(4) = xccorn(1) + dxc
         yccorn(4) = yccorn(1)
         call mapc2p(xccorn(4),yccorn(4),xpcorn(4),ypcorn(4))

         ! compute normals to left and bottom side:
         ax =  (ypcorn(2) - ypcorn(1))
         ay = -(xpcorn(2) - xpcorn(1))
         anorm = sqrt(ax*ax + ay*ay)
         aux(1,i,j) = ax/anorm
         aux(2,i,j) = ay/anorm
         aux(3,i,j) = anorm/dyc

         bx = -(ypcorn(4) - ypcorn(1))
         by =  (xpcorn(4) - xpcorn(1))
         bnorm = sqrt(bx*bx + by*by)
         aux(4,i,j) = bx/bnorm
         aux(5,i,j) = by/bnorm
         aux(6,i,j) = bnorm/dxc

         ! compute area of physical cell from four corners:
         xpcorn(5) = xpcorn(1)
         ypcorn(5) = ypcorn(1)

         area = 0.0d0
         do ic=1,4
            area = area + 0.5d0 * (ypcorn(ic)+ypcorn(ic+1)) * &
                                  (xpcorn(ic+1)-xpcorn(ic))
         enddo
         aux(7,i,j) = area / (dxc*dyc)

        enddo
    enddo

    return
end subroutine setaux

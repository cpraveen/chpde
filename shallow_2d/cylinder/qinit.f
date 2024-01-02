!     =====================================================
      SUBROUTINE qinit (meqn, mbc, mx, my, xlower, ylower, dx, dy, q,
     &                  maux, aux)
!     =====================================================
! 
!     # Set initial conditions for q.
! 
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      COMMON /comic/ hl,hr,hul,hur
! 
      DO i=1,mx
         xc=xlower+(i-0.5d0)*dx
         xclow=xlower+(i-1.0d0)*dx
         DO j=1,my
            yc=ylower+(j-0.5d0)*dy
            yclow=ylower+(j-1.0d0)*dy
  
!           map the center of this computational cell to physical
!           coordinates before evaluating the initial value funcion:
            CALL mapc2p (xc, yc, xp, yp)
  
!          q(i,j,1) = 1.d0  +
!    &                3.d0*dexp(-50.d0*((xp-1.0d0)**2 + (yp-1.0d0)**2))
!    &                dexp(-100.d0*((xp-0.8d0)**2 + (yp-1.0d0)**2))
  
            CALL cellave (xclow, yclow, dx, dy, win)
            q(1,i,j)=hl*win+hr*(1.d0-win)
            q(2,i,j)=hul*win+hur*(1.d0-win)
            q(3,i,j)=0.0d0
  
         END DO
      END DO
  
      RETURN
      END

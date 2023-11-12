! =========================================================
      subroutine qinit (meqn, mbc, mx, xlower, dx, q, maux, aux)
! =========================================================
! 
!     # Set initial conditions for q.
! 
! 
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
! 
! 
      do i=1,mx
         xcell=xlower+(i-0.5d0)*dx
!        q(1,i) = g0((xcell-50.d0)/5.d0)
         q(1,i)=0.d0
         q(2,i)=0.d0
      end do
! 
      return
      end

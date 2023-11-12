!     ============================================
      subroutine setaux (mbc, mx, xlower, dx, maux, aux)
!     ============================================
! 
!     # set auxiliary arrays
!     # variable coefficient acoustics
!     #  aux(1,i) = impedance Z in i'th cell
!     #  aux(2,i) = sound speed c in i'th cell
! 
! 
!     # rho value for i-1 < x < i is determined by rho(i)
!     # input from setprob.rho in setprob.f
! 
! 
      implicit double precision (a-h,o-z)
      dimension aux(maux,1-mbc:mx+mbc)
      common /comaux/ rho(1000)
  
  
      do i=1-mbc,mx+mbc
         xcell=xlower+(i-0.5d0)*dx
!        # truncate to an integer to determine which layer xcell is in:
         ix=xcell+1
         if (ix.lt.1) ix=1
         aux(2,i)=1.d0
         aux(1,i)=rho(ix)*aux(2,i)
!        # uniform for comparison:
!        aux(1,i) = aux(2,i)
      end do
  
      open (unit=31,file='fort.aux',status='unknown',form='formatted')
  
      do i=1,mx
         write (31,10) aux(1,i),aux(2,i)
   10    format (2e16.6)
      end do
  
      close (unit=31)
! 
      return
      end

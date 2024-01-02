!     =====================================================
      SUBROUTINE rpt2 (ixy, imp, maxm, meqn, mwaves, maux, mbc, mx, 
     &                 ql, qr, aux1, aux2, aux3, asdq, bmasdq, bpasdq)
!     =====================================================
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
! 
!     # Riemann solver in the transverse direction for the shallow water
!     # equations  on a quadrilateral grid.
! 
!     # Split asdq (= A^* \Delta q, where * = + or -)
!     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
!     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
! 
!     # Uses Roe averages and other quantities which were
!     # computed in rpn2sh and stored in the common block comroe.
! 
      DIMENSION ql(meqn,1-mbc:maxm+mbc)
      DIMENSION qr(meqn,1-mbc:maxm+mbc)
      DIMENSION asdq(meqn,1-mbc:maxm+mbc)
      DIMENSION bmasdq(meqn,1-mbc:maxm+mbc)
      DIMENSION bpasdq(meqn,1-mbc:maxm+mbc)
      DIMENSION aux1(maux,1-mbc:maxm+mbc)
      DIMENSION aux2(maux,1-mbc:maxm+mbc)
      DIMENSION aux3(maux,1-mbc:maxm+mbc)
! 
      DIMENSION u(-1:maxm+2), v(-1:maxm+2), a(-1:maxm+2), h(-1:maxm+2)
      DIMENSION alf(-1:maxm+2)
      DIMENSION beta(-1:maxm+2)
      DIMENSION wave(3,3,-1:maxm+2)
      DIMENSION s(3,-1:maxm+2)
      DIMENSION delta(3)
      COMMON /sw/ g
! 
      IF (ixy.eq.1) THEN
         inx=4
         iny=5
         ilenrat=6
      ELSE
         inx=1
         iny=2
         ilenrat=3
      END IF
! 
!     imp is used to flag whether wave is going to left or right,
!     since states and grid are different on each side
! 
      IF (imp.eq.1) THEN
!        asdq = amdq, moving to left
         ix1=2-mbc
         ixm1=mx+mbc
      ELSE
!        asdq = apdq, moving to right
         ix1=1-mbc
         ixm1=mx+mbc
      END IF
! 
!        --------------
!        # up-going:
!        --------------
! 
  
!       # determine rotation matrix for interface above cell, using aux3
!               [ alf  beta ]
!               [-beta  alf ]
! 
      DO i=ix1,ixm1
! 
         IF (imp.eq.1) THEN
            i1=i-1
         ELSE
            i1=i
         END IF
! 
         alf(i)=aux3(inx,i1)
         beta(i)=aux3(iny,i1)
         h(i)=ql(1,i1)
         u(i)=(alf(i)*ql(2,i1)+beta(i)*ql(3,i1))/h(i)
         v(i)=(-beta(i)*ql(2,i1)+alf(i)*ql(3,i1))/h(i)
         a(i)=dsqrt(g*h(i))
      END DO
! 
! 
!     # now split asdq into waves:
! 
!     # find a1 thru a3, the coefficients of the 3 eigenvectors:
      DO 10 i=ix1,ixm1
         delta(1)=asdq(1,i)
         delta(2)=alf(i)*asdq(2,i)+beta(i)*asdq(3,i)
         delta(3)=-beta(i)*asdq(2,i)+alf(i)*asdq(3,i)
         a1=((u(i)+a(i))*delta(1)-delta(2))*(0.50d0/a(i))
         a2=-v(i)*delta(1)+delta(3)
         a3=(-(u(i)-a(i))*delta(1)+delta(2))*(0.50d0/a(i))
! 
!        # Compute the waves.
! 
         wave(1,1,i)=a1
         wave(2,1,i)=alf(i)*a1*(u(i)-a(i))-beta(i)*a1*v(i)
         wave(3,1,i)=beta(i)*a1*(u(i)-a(i))+alf(i)*a1*v(i)
         s(1,i)=(u(i)-a(i))*aux3(ilenrat,i1)
! 
         wave(1,2,i)=0.0d0
         wave(2,2,i)=-beta(i)*a2
         wave(3,2,i)=alf(i)*a2
         s(2,i)=u(i)*aux3(ilenrat,i1)
! 
         wave(1,3,i)=a3
         wave(2,3,i)=alf(i)*a3*(u(i)+a(i))-beta(i)*a3*v(i)
         wave(3,3,i)=beta(i)*a3*(u(i)+a(i))+alf(i)*a3*v(i)
         s(3,i)=(u(i)+a(i))*aux3(ilenrat,i1)
   10 CONTINUE
! 
! 
!    # compute flux difference bpasdq
!    --------------------------------
! 
      DO 30 m=1,3
         DO 30 i=ix1,ixm1
         bpasdq(m,i)=0.d0
         DO 20 mw=1,mwaves
            bpasdq(m,i)=bpasdq(m,i)+dmax1(s(mw,i),0.d0)*wave(m,mw,i)
   20    CONTINUE
   30 CONTINUE
! 
! 
!        --------------
!        # down-going:
!        --------------
! 
  
!       # determine rotation matrix for interface below cell, using aux2
!               [ alf  beta ]
!               [-beta  alf ]
! 
      DO i=ix1,ixm1
! 
         IF (imp.eq.1) THEN
            i1=i-1
         ELSE
            i1=i
         END IF
! 
         alf(i)=aux2(inx,i1)
         beta(i)=aux2(iny,i1)
         u(i)=(alf(i)*ql(2,i1)+beta(i)*ql(3,i1))/h(i)
         v(i)=(-beta(i)*ql(2,i1)+alf(i)*ql(3,i1))/h(i)
      END DO
! 
! 
!     # now split asdq into waves:
! 
!     # find a1 thru a3, the coefficients of the 3 eigenvectors:
      DO 40 i=ix1,ixm1
         delta(1)=asdq(1,i)
         delta(2)=alf(i)*asdq(2,i)+beta(i)*asdq(3,i)
         delta(3)=-beta(i)*asdq(2,i)+alf(i)*asdq(3,i)
         a1=((u(i)+a(i))*delta(1)-delta(2))*(0.50d0/a(i))
         a2=-v(i)*delta(1)+delta(3)
         a3=(-(u(i)-a(i))*delta(1)+delta(2))*(0.50d0/a(i))
! 
!        # Compute the waves.
! 
         wave(1,1,i)=a1
         wave(2,1,i)=alf(i)*a1*(u(i)-a(i))-beta(i)*a1*v(i)
         wave(3,1,i)=beta(i)*a1*(u(i)-a(i))+alf(i)*a1*v(i)
         s(1,i)=(u(i)-a(i))*aux2(ilenrat,i1)
! 
         wave(1,2,i)=0.0d0
         wave(2,2,i)=-beta(i)*a2
         wave(3,2,i)=alf(i)*a2
         s(2,i)=u(i)*aux2(ilenrat,i1)
! 
         wave(1,3,i)=a3
         wave(2,3,i)=alf(i)*a3*(u(i)+a(i))-beta(i)*a3*v(i)
         wave(3,3,i)=beta(i)*a3*(u(i)+a(i))+alf(i)*a3*v(i)
         s(3,i)=(u(i)+a(i))*aux2(ilenrat,i1)
   40 CONTINUE
! 
! 
!    # compute flux difference bmasdq
!    --------------------------------
! 
      DO 60 m=1,3
         DO 60 i=ix1,ixm1
         bmasdq(m,i)=0.d0
         DO 50 mw=1,mwaves
            bmasdq(m,i)=bmasdq(m,i)+dmin1(s(mw,i),0.d0)*wave(m,mw,i)
   50    CONTINUE
   60 CONTINUE
! 
! 
      RETURN
      END

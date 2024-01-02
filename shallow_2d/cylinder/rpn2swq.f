  
! 
! 
!     =====================================================
      SUBROUTINE rpn2 (ixy, maxm, meqn, mwaves, maux, mbc, mx, ql, qr, 
     &                 auxl, auxr, wave, s, amdq, apdq)
!     =====================================================
! 
!     # Roe-solver for the 2D shallow water equations
!     #  on a quadrilateral grid
! 
!     # solve Riemann problems along one slice of data.
! 
!     # On input, ql contains the state vector at the left edge of each
!     #           qr contains the state vector at the right edge of each
! 
!     # This data is along a slice in the x-direction if ixy=1
!     #                            or the y-direction if ixy=2.
!     # On output, wave contains the waves, s the speeds,
!     # and amdq, apdq the decomposition of the flux difference
!     #   f(qr(i-1)) - f(ql(i))
!     # into leftgoing and rightgoing parts respectively.
!     # With the Roe solver we have
!     #    amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
!     # where A is the Roe matrix.  An entropy fix can also be incorpora
!     # into the flux differences.
! 
!     # Note that the i'th Riemann problem has left state qr(i-1,:)
!     #                                    and right state ql(i,:)
!     # From the basic clawpack routines, this routine is called with ql
! 
! 
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
! 
      DIMENSION wave(meqn,mwaves,1-mbc:maxm+mbc)
      DIMENSION s(mwaves,1-mbc:maxm+mbc)
      DIMENSION ql(meqn,1-mbc:maxm+mbc)
      DIMENSION qr(meqn,1-mbc:maxm+mbc)
      DIMENSION apdq(meqn,1-mbc:maxm+mbc)
      DIMENSION amdq(meqn,1-mbc:maxm+mbc)
      DIMENSION auxl(maux,1-mbc:maxm+mbc)
      DIMENSION auxr(maux,1-mbc:maxm+mbc)
! 
!     local arrays
!     ------------
      DIMENSION delta(3)
      LOGICAL efix
      DIMENSION unorl(-1:maxm+2), unorr(-1:maxm+2)
      DIMENSION utanl(-1:maxm+2), utanr(-1:maxm+2)
      DIMENSION alf(-1:maxm+2)
      DIMENSION beta(-1:maxm+2)
      DIMENSION u(-1:maxm+2)
      DIMENSION v(-1:maxm+2)
      DIMENSION a(-1:maxm+2)
      DIMENSION h(-1:maxm+2)
  
      COMMON /sw/ g
! 
      DATA efix /.true./
! 
! 
!     # rotate the velocities q(2) and q(3) so that it is aligned with g
!     # normal.  The normal vector for the face at the i'th Riemann prob
!     # is stored in the aux array
!     # in locations (1,2) if ixy=1 or (4,5) if ixy=2.  The ratio of the
!     # length of the cell side to the length of the computational cell
!     # is stored in aux(3) or aux(6) respectively.
! 
! 
      IF (ixy.eq.1) THEN
         inx=1
         iny=2
         ilenrat=3
      ELSE
         inx=4
         iny=5
         ilenrat=6
      END IF
! 
!       # determine rotation matrix
!               [ alf  beta ]
!               [-beta  alf ]
! 
!       # note that this reduces to identity on standard cartesian grid
! 
      DO i=2-mbc,mx+mbc
         alf(i)=auxl(inx,i)
         beta(i)=auxl(iny,i)
         unorl(i)=alf(i)*ql(2,i)+beta(i)*ql(3,i)
         unorr(i-1)=alf(i)*qr(2,i-1)+beta(i)*qr(3,i-1)
         utanl(i)=-beta(i)*ql(2,i)+alf(i)*ql(3,i)
         utanr(i-1)=-beta(i)*qr(2,i-1)+alf(i)*qr(3,i-1)
      END DO
! 
! 
!     # compute the Roe-averaged variables needed in the Roe solver.
!     # These are stored in the common block comroe since they are
!     # later used in routine rpt2 to do the transverse wave splitting.
! 
      DO 10 i=2-mbc,mx+mbc
         h(i)=(qr(1,i-1)+ql(1,i))*0.50d0
         hsqrtl=dsqrt(qr(1,i-1))
         hsqrtr=dsqrt(ql(1,i))
         hsq2=hsqrtl+hsqrtr
         u(i)=(unorr(i-1)/hsqrtl+unorl(i)/hsqrtr)/hsq2
         v(i)=(utanr(i-1)/hsqrtl+utanl(i)/hsqrtr)/hsq2
         a(i)=dsqrt(g*h(i))
   10 CONTINUE
! 
! 
!     # now split the jump in q at each interface into waves
! 
!     # find a1 thru a3, the coefficients of the 3 eigenvectors:
      DO 20 i=2-mbc,mx+mbc
         delta(1)=ql(1,i)-qr(1,i-1)
         delta(2)=unorl(i)-unorr(i-1)
         delta(3)=utanl(i)-utanr(i-1)
         a1=((u(i)+a(i))*delta(1)-delta(2))*(0.50d0/a(i))
         a2=-v(i)*delta(1)+delta(3)
         a3=(-(u(i)-a(i))*delta(1)+delta(2))*(0.50d0/a(i))
! 
!        # Compute the waves.
! 
         wave(1,1,i)=a1
         wave(1,2,i)=alf(i)*a1*(u(i)-a(i))-beta(i)*a1*v(i)
         wave(1,3,i)=beta(i)*a1*(u(i)-a(i))+alf(i)*a1*v(i)
         s(1,i)=(u(i)-a(i))*auxl(ilenrat,i)
! 
         wave(2,1,i)=0.0d0
         wave(2,2,i)=-beta(i)*a2
         wave(2,3,i)=alf(i)*a2
         s(2,i)=u(i)*auxl(ilenrat,i)
! 
         wave(3,1,i)=a3
         wave(3,2,i)=alf(i)*a3*(u(i)+a(i))-beta(i)*a3*v(i)
         wave(3,3,i)=beta(i)*a3*(u(i)+a(i))+alf(i)*a3*v(i)
         s(3,i)=(u(i)+a(i))*auxl(ilenrat,i)
   20 CONTINUE
! 
! 
!    # compute flux differences amdq and apdq.
!    ---------------------------------------
! 
      IF (efix) GO TO 50
! 
!     # no entropy fix
!     ----------------
! 
!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves
! 
      DO 40 m=1,3
         DO 40 i=2-mbc,mx+mbc
         amdq(m,i)=0.d0
         apdq(m,i)=0.d0
         DO 30 mw=1,mwaves
            IF (s(mw,i).lt.0.d0) THEN
               amdq(m,i)=amdq(m,i)+s(mw,i)*wave(m,mw,i)
            ELSE
               apdq(m,i)=apdq(m,i)+s(mw,i)*wave(m,mw,i)
            END IF
   30    CONTINUE
   40 CONTINUE
      GO TO 130
! 
!-----------------------------------------------------
! 
   50 CONTINUE
! 
!     # With entropy fix
!     ------------------
! 
!    # compute flux differences amdq and apdq.
!    # First compute amdq as sum of s*wave for left going waves.
!    # Incorporate entropy fix by adding a modified fraction of wave
!    # if s should change sign.
! 
      DO 100 i=2-mbc,mx+mbc
!        check 1-wave
         him1=qr(1,i-1)
         s0=(unorr(i-1)/him1-dsqrt(g*him1))*auxl(ilenrat,i)
!        check for fully supersonic case :
         IF (s0.gt.0.0d0.and.s(1,i).gt.0.0d0) THEN
            DO 60 m=1,3
               amdq(m,i)=0.0d0
   60       CONTINUE
            GO TO 100
         END IF
! 
         h1=qr(1,i-1)+wave(1,1,i)
         hu1=unorr(i-1)+alf(i)*wave(2,1,i)+beta(i)*wave(3,1,i)
         s1=(hu1/h1-dsqrt(g*h1))*auxl(ilenrat,i)
         IF (s0.lt.0.0d0.and.s1.gt.0.0d0) THEN
!           transonic rarefaction in 1-wave
            sfract=s0*((s1-s(1,i))/(s1-s0))
         ELSE IF (s(1,i).lt.0.0d0) THEN
!           1-wave is leftgoing
            sfract=s(1,i)
         ELSE
!           1-wave is rightgoing
            sfract=0.0d0
         END IF
         DO 70 m=1,3
            amdq(m,i)=sfract*wave(m,1,i)
   70    CONTINUE
!           check 2-wave
         IF (s(2,i).gt.0.0d0) THEN
!           2 and 3 waves are right-going
            GO TO 100
         END IF
  
         DO 80 m=1,3
            amdq(m,i)=amdq(m,i)+s(2,i)*wave(m,2,i)
   80    CONTINUE
! 
!        check 3-wave
! 
         hi=ql(1,i)
         s03=(unorl(i)/hi+dsqrt(g*hi))*auxl(ilenrat,i)
         h3=ql(1,i)-wave(1,3,i)
         hu3=unorl(i)-(alf(i)*wave(2,3,i)+beta(i)*wave(3,3,i))
         s3=(hu3/h3+dsqrt(g*h3))*auxl(ilenrat,i)
         IF (s3.lt.0.0d0.and.s03.gt.0.0d0) THEN
!           transonic rarefaction in 3-wave
            sfract=s3*((s03-s(3,i))/(s03-s3))
         ELSE IF (s(3,i).lt.0.0d0) THEN
!           3-wave is leftgoing
            sfract=s(3,i)
         ELSE
!           3-wave is rightgoing
            GO TO 100
         END IF
         DO 90 m=1,3
            amdq(m,i)=amdq(m,i)+sfract*wave(m,3,i)
   90    CONTINUE
  100 CONTINUE
! 
!     compute rightgoing flux differences :
! 
      DO 120 m=1,3
         DO 120 i=2-mbc,mx+mbc
         df=0.0d0
         DO 110 mw=1,mwaves
            df=df+s(mw,i)*wave(m,mw,i)
  110    CONTINUE
         apdq(m,i)=df-amdq(m,i)
  120 CONTINUE
! 
! 
  130 CONTINUE
      RETURN
      END

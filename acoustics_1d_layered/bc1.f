  
!     =================================================================
      subroutine bc1 (meqn, mbc, mx, xlower, dx, q, maux, aux, t, dt,
     &mthbc)
!     =================================================================
! 
!     # Standard boundary condition choices for claw2
! 
!     # Modified for layered medium:
!     #   For small t oscillating solid wall BCs are used to generate pu
!     #   For t>50, switch to periodic boundary conditions.  The code
!     #     can then be run to much larger t to observe how the pulse
!     #     behaves as it loops around.
! 
!     # At each boundary  k = 1 (left),  2 (right):
!     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
!     #            =  1  for zero-order extrapolation
!     #            =  2  for periodic boundary coniditions
!     #            =  3  for solid walls, assuming this can be implement
!     #                  by reflecting the data about the boundary and t
!     #                  negating the 2'nd component of q.
!     ------------------------------------------------
! 
!     # Extend the data from the computational region
!     #      i = 1, 2, ..., mx2
!     # to the virtual cells outside the region, with
!     #      i = 1-ibc  and   i = mx+ibc   for ibc=1,...,mbc
! 
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
  
      dimension mthbc(2)
      common /comwall/ pi,t1,a1,tw1,t2,a2,tw2
  
!     # switch to periodic boundary conditions once pulse is generated:
      if (t.gt.50.d0.and.t.lt.51.d0) then
         mthbc(1)=2
         mthbc(2)=2
         do m=1,meqn
            do ibc=1,mbc
               aux(m,1-ibc)=aux(m,mx+1-ibc)
            end do
         end do
      end if
  
!-------------------------------------------------------
!     # left boundary:
!-------------------------------------------------------
      go to (10,40,60,80),mthbc(1)+1
! 
   10 continue
!     # user-specified boundary conditions
!     # oscillating wall
      do 20 m=1,meqn
         do 20 ibc=1,mbc
         q(m,1-ibc)=q(m,ibc)
   20 continue
!     # wall velocity:
      vwall=a1*g0((t-t1)/tw1)
!     # adjust the normal velocity:
      do 30 ibc=1,mbc
         q(2,1-ibc)=2.0d0*vwall-q(2,ibc)
   30 continue
      go to 110
! 
   40 continue
!     # zero-order extrapolation:
      do 50 m=1,meqn
         do 50 ibc=1,mbc
         q(m,1-ibc)=q(m,1)
   50 continue
      go to 110
  
   60 continue
!     # periodic:
      do 70 m=1,meqn
         do 70 ibc=1,mbc
         q(m,1-ibc)=q(m,mx+1-ibc)
   70 continue
      go to 110
  
   80 continue
!     # solid wall (assumes 2'nd component is velocity or momentum in x)
      do 90 m=1,meqn
         do 90 ibc=1,mbc
         q(m,1-ibc)=q(m,ibc)
   90 continue
!     # negate the normal velocity:
      do 100 ibc=1,mbc
         q(2,1-ibc)=-q(2,ibc)
  100 continue
      go to 110
  
  110 continue
  
!-------------------------------------------------------
!     # right boundary:
!-------------------------------------------------------
      go to (120,150,170,190),mthbc(2)+1
! 
  120 continue
!     # user-specified boundary conditions
!     # oscillating wall
      do 130 m=1,meqn
         do 130 ibc=1,mbc
         q(m,mx+ibc)=q(m,mx+1-ibc)
  130 continue
!     # wall velocity:
      vwall=a2*g0((t-t2)/tw2)
!     # adjust the normal velocity:
      do 140 ibc=1,mbc
         q(2,mx+ibc)=2.0d0*vwall-q(2,mx+1-ibc)
  140 continue
      go to 220
  
  150 continue
!     # zero-order extrapolation:
      do 160 m=1,meqn
         do 160 ibc=1,mbc
         q(m,mx+ibc)=q(m,mx)
  160 continue
      go to 220
  
  170 continue
!     # periodic:
      do 180 m=1,meqn
         do 180 ibc=1,mbc
         q(m,mx+ibc)=q(m,ibc)
  180 continue
      go to 220
  
  190 continue
!     # solid wall (assumes 2'nd component is velocity or momentum in x)
      do 200 m=1,meqn
         do 200 ibc=1,mbc
         q(m,mx+ibc)=q(m,mx+1-ibc)
  200 continue
      do 210 ibc=1,mbc
         q(2,mx+ibc)=-q(2,mx+1-ibc)
  210 continue
      go to 220
  
  220 continue
! 
      return
      end
  
  
!     ===============================
      double precision function g0 (t)
!     ===============================
  
      implicit double precision (a-h,o-z)
      common /comwall/ pi,t1,a1,t2,a2,tw1,tw2
  
      if (dabs(t).lt.1.d0) then
         g0=1.d0+dcos(pi*t)
      else
         g0=0.d0
      end if
  
      return
      end

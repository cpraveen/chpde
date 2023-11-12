  
!     =================================================================
      subroutine bc1 (meqn, mbc, mx, xlower, dx, q, maux, aux, t, dt,
     &                mthbc)
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
      go to (10,20,30,40) mthbc(1)+1
! 
   10 continue
!     # user-specified boundary conditions
!     # oscillating wall
      do m=1,meqn
         do ibc=1,mbc
            q(m,1-ibc)=q(m,ibc)
         end do
      end do
!     # wall velocity:
      vwall=a1*g0((t-t1)/tw1)
!     # adjust the normal velocity:
      do ibc=1,mbc
         q(2,1-ibc)=2.0d0*vwall-q(2,ibc)
      end do
      go to 50
! 
   20 continue
!     # zero-order extrapolation:
      do m=1,meqn
         do ibc=1,mbc
            q(m,1-ibc)=q(m,1)
         end do
      end do
      go to 50
  
   30 continue
!     # periodic:
      do m=1,meqn
         do ibc=1,mbc
            q(m,1-ibc)=q(m,mx+1-ibc)
         end do
      end do
      go to 50
  
   40 continue
!     # solid wall (assumes 2'nd component is velocity or momentum in x)
      do m=1,meqn
         do ibc=1,mbc
            q(m,1-ibc)=q(m,ibc)
         end do
      end do
!     # negate the normal velocity:
      do ibc=1,mbc
         q(2,1-ibc)=-q(2,ibc)
      end do
      go to 50
  
   50 continue
  
!-------------------------------------------------------
!     # right boundary:
!-------------------------------------------------------
      go to (60,70,80,90) mthbc(2)+1
! 
   60 continue
!     # user-specified boundary conditions
!     # oscillating wall
      do m=1,meqn
         do ibc=1,mbc
            q(m,mx+ibc)=q(m,mx+1-ibc)
         end do
      end do
!     # wall velocity:
      vwall=a2*g0((t-t2)/tw2)
!     # adjust the normal velocity:
      do ibc=1,mbc
         q(2,mx+ibc)=2.0d0*vwall-q(2,mx+1-ibc)
      end do
      go to 100
  
   70 continue
!     # zero-order extrapolation:
      do m=1,meqn
         do ibc=1,mbc
            q(m,mx+ibc)=q(m,mx)
         end do
      end do
      go to 100
  
   80 continue
!     # periodic:
      do m=1,meqn
         do ibc=1,mbc
            q(m,mx+ibc)=q(m,ibc)
         end do
      end do
      go to 100
  
   90 continue
!     # solid wall (assumes 2'nd component is velocity or momentum in x)
      do m=1,meqn
         do ibc=1,mbc
            q(m,mx+ibc)=q(m,mx+1-ibc)
         end do
      end do
      do ibc=1,mbc
         q(2,mx+ibc)=-q(2,mx+1-ibc)
      end do
      go to 100
  
  100 continue
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

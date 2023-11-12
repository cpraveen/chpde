! ===================================================================
      subroutine step1 (meqn, mwaves, mbc, maux, mx, q, aux, dx, dt,
     &method, mthlim, cfl, f, wave, s, amdq, apdq, dtdx, use_fwave, rp1)
! ===================================================================
! 
!     # Take one time step, updating q.
! 
!     method(1) = 1   ==>  Godunov method
!     method(1) = 2   ==>  Slope limiter method
!     mthlim(p)  controls what limiter is used in the pth family
! 
! 
!     amdq, apdq, wave, s, and f are used locally:
! 
!     amdq(1-mbc:mx+mbc, meqn) = left-going flux-differences
!     apdq(1-mbc:mx+mbc, meqn) = right-going flux-differences
!        e.g. amdq(i,m) = m'th component of A^- \Delta q from i'th Riema
!                         problem (between cells i-1 and i).
! 
!     wave(1-mbc:mx+mbc, meqn, mwaves) = waves from solution of
!                                           Riemann problems,
!            wave(i,m,mw) = mth component of jump in q across
!                           wave in family mw in Riemann problem between
!                           states i-1 and i.
! 
!     s(1-mbc:mx+mbc, mwaves) = wave speeds,
!            s(i,mw) = speed of wave in family mw in Riemann problem bet
!                      states i-1 and i.
! 
!     f(1-mbc:mx+mbc, meqn) = correction fluxes for second order method
!            f(i,m) = mth component of flux at left edge of ith cell
!     ------------------------------------------------------------------
! 
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
      dimension f(meqn,1-mbc:mx+mbc)
      dimension s(mwaves,1-mbc:mx+mbc)
      dimension wave(meqn,mwaves,1-mbc:mx+mbc)
      dimension amdq(meqn,1-mbc:mx+mbc)
      dimension apdq(meqn,1-mbc:mx+mbc)
      dimension dtdx(1-mbc:mx+mbc)
      dimension method(7), mthlim(mwaves)
      logical use_fwave
      logical limit
      external rp1
      common /comlim/ mylim
! 
!     # check if any limiters are used:
      limit=.false.
      do mw=1,mwaves
         if (mthlim(mw).gt.0) limit=.true.
      end do
! 
      mcapa=method(6)
      do i=1-mbc,mx+mbc
         if (mcapa.gt.0) then
            if (aux(i,mcapa).le.0.d0) then
               write (6,*) 'Error -- capa must be positive'
               stop
            end if
            dtdx(i)=dt/(dx*aux(mcapa,i))
         else
            dtdx(i)=dt/dx
         end if
      end do
! 
! 
!     # solve Riemann problem at each interface
!     -----------------------------------------
! 
      call rp1 (mx, meqn, mwaves, maux, mbc, mx, q, q, aux, aux, wave,
     &s, amdq, apdq)
! 
!     # Modify q for Godunov update:
!     # Note this may not correspond to a conservative flux-differencing
!     # for equations not in conservation form.  It is conservative if
!     # amdq + apdq = f(q(i)) - f(q(i-1)).
! 
      do i=1,mx+1
         do m=1,meqn
            q(m,i)=q(m,i)-dtdx(i)*apdq(m,i)
            q(m,i-1)=q(m,i-1)-dtdx(i-1)*amdq(m,i)
         end do
      end do
  
! 
!     # compute maximum wave speed:
      cfl=0.d0
      do mw=1,mwaves
         do i=1,mx+1
!          # if s>0 use dtdx(i) to compute CFL,
!          # if s<0 use dtdx(i-1) to compute CFL:
            cfl=dmax1(cfl,dtdx(i)*s(mw,i),-dtdx(i-1)*s(mw,i))
         end do
      end do
! 
      if (method(2).eq.1) go to 30
! 
!     # compute correction fluxes for second order q_{xx} terms:
!     ----------------------------------------------------------
! 
      do m=1,meqn
         do i=1-mbc,mx+mbc
            f(m,i)=0.d0
         end do
      end do
! 
!     # apply limiter to waves:
      if (limit) then
         if (mylim.eq.1) then
!          # transmission-based limiter:
            call trlimit (meqn, mwaves, maux, mbc, mx, aux, wave, s,
     &       mthlim)
         else
!          # original:
            call limiter (mx, meqn, mwaves, mbc, mx, wave, s, mthlim)
         end if
      end if
  
! 
      do 20 i=1,mx+1
         do 20 m=1,meqn
         do 10 mw=1,mwaves
            dtdxave=0.5d0*(dtdx(i-1)+dtdx(i))
            f(m,i)=f(m,i)+0.5d0*dabs(s(mw,i))*(1.d0-dabs(s(mw,i))*
     &       dtdxave)*wave(m,mw,i)
! 
!              # third order corrections:
!              # (still experimental... works well for smooth solutions
!              # with no limiters but not well with limiters so far.
! 
  
            if (method(2).lt.3) go to 10
            if (s(i,mw).gt.0.d0) then
               dq2=wave(m,mw,i)-wave(m,mw,i-1)
            else
               dq2=wave(m,mw,i+1)-wave(m,mw,i)
            end if
            f(m,i)=f(m,i)-s(mw,i)/6.d0*(1.d0-(s(mw,i)*dtdxave)**2)*dq2
  
   10    continue
   20 continue
! 
! 
!     # update q by differencing correction fluxes
!     ============================================
! 
!     # (Note:  Godunov update has already been performed above)
! 
      do m=1,meqn
         do i=1,mx
            q(m,i)=q(m,i)-dtdx(i)*(f(m,i+1)-f(m,i))
         end do
      end do
! 
   30 continue
      return
      end

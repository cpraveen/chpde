c
c
c ===================================================================
      subroutine step1(meqn,mwaves,mbc,maux,mx,q,aux,dx,dt,method,
     &                 mthlim,cfl,f,wave,s,amdq,apdq,dtdx,use_fwave,rp1)
c ===================================================================
c
c     # Take one time step, updating q.
c
c     method(1) = 1   ==>  Godunov method
c     method(1) = 2   ==>  Slope limiter method
c     mthlim(p)  controls what limiter is used in the pth family
c
c
c     amdq, apdq, wave, s, and f are used locally:
c
c     amdq(1-mbc:mx+mbc, meqn) = left-going flux-differences
c     apdq(1-mbc:mx+mbc, meqn) = right-going flux-differences
c        e.g. amdq(i,m) = m'th component of A^- \Delta q from i'th Riemann
c                         problem (between cells i-1 and i).
c
c     wave(1-mbc:mx+mbc, meqn, mwaves) = waves from solution of
c                                           Riemann problems,
c            wave(i,m,mw) = mth component of jump in q across
c                           wave in family mw in Riemann problem between
c                           states i-1 and i.
c
c     s(1-mbc:mx+mbc, mwaves) = wave speeds,
c            s(i,mw) = speed of wave in family mw in Riemann problem between
c                      states i-1 and i.
c
c     f(1-mbc:mx+mbc, meqn) = correction fluxes for second order method
c            f(i,m) = mth component of flux at left edge of ith cell 
c     --------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      dimension    q(meqn,1-mbc:mx+mbc)
      dimension  aux(maux,1-mbc:mx+mbc)
      dimension    f(meqn,1-mbc:mx+mbc)
      dimension    s(mwaves,1-mbc:mx+mbc)
      dimension wave(meqn,mwaves,1-mbc:mx+mbc)
      dimension amdq(meqn,1-mbc:mx+mbc)
      dimension apdq(meqn,1-mbc:mx+mbc)
      dimension dtdx(1-mbc:mx+mbc)
      dimension method(7),mthlim(mwaves)
      logical use_fwave
      logical limit
      external rp1
      common /comlim/ mylim,mrplim(2)
c
c     # check if any limiters are used:
      limit = .false.
      do mw=1,mwaves
         if (mylim.eq.0 .and. mthlim(mw) .gt. 0) limit = .true.
         if (mylim.eq.1 .and. mrplim(mw) .gt. 0) limit = .true.
      enddo
c
      mcapa = method(6)
      do i=1-mbc,mx+mbc
         if (mcapa.gt.0) then
           if (aux(i,mcapa) .le. 0.d0) then
		write(6,*) 'Error -- capa must be positive'
		stop
           endif
             dtdx(i) = dt / (dx*aux(mcapa,i))
	    else
             dtdx(i) = dt/dx
        endif
      enddo
c
c
c
c     # solve Riemann problem at each interface 
c     -----------------------------------------
c
      call rp1(mx,meqn,mwaves,maux,mbc,mx,q,q,aux,aux,wave,s,amdq,apdq)
c
c     # Modify q for Godunov update:
c     # Note this may not correspond to a conservative flux-differencing
c     # for equations not in conservation form.  It is conservative if
c     # amdq + apdq = f(q(i)) - f(q(i-1)).
c
      do i=1,mx+1
         do m=1,meqn
            q(m,i) = q(m,i) - dtdx(i)*apdq(m,i)
            q(m,i-1) = q(m,i-1) - dtdx(i-1)*amdq(m,i)
         enddo
      enddo

c
c     # compute maximum wave speed:
      cfl = 0.d0
      do mw=1,mwaves
	 do i=1,mx+1
c          # if s>0 use dtdx(i) to compute CFL,
c          # if s<0 use dtdx(i-1) to compute CFL:
	   cfl = dmax1(cfl, dtdx(i)*s(mw,i), -dtdx(i-1)*s(mw,i))
       enddo
      enddo
c
      if (method(2) .eq. 1) go to 900
c
c     # compute correction fluxes for second order q_{xx} terms:
c     ----------------------------------------------------------
c
      do m = 1, meqn
            do i = 1-mbc, mx+mbc
               f(m,i) = 0.d0
            enddo
      enddo
c
c     # apply limiter to waves:
      if (limit) then
        if (mylim .eq. 1) then
c          # transmission-based limiter:
           call trlimit(meqn,mwaves,maux,mbc,mx,aux,wave,s,mrplim)
         else
c          # original:
           call limiter(mx,meqn,mwaves,mbc,mx,wave,s,mthlim)
         endif
        endif

c
      do 120 i=1,mx+1
        do 120 m=1,meqn
	    do 110 mw=1,mwaves
	       dtdxave = 0.5d0 * (dtdx(i-1) + dtdx(i))
	       f(m,i) = f(m,i) + 0.5d0 * dabs(s(mw,i))
     &                * (1.d0 - dabs(s(mw,i))*dtdxave) * wave(m,mw,i)
c
c              # third order corrections:
c              # (still experimental... works well for smooth solutions
c              # with no limiters but not well with limiters so far.
c

               if (method(2).lt.3) go to 110
               if (s(i,mw) .gt. 0.d0) then
                   dq2 = wave(m,mw,i) - wave(m,mw,i-1)
                 else
                   dq2 = wave(m,mw,i+1) - wave(m,mw,i)
                 endif
               f(m,i) = f(m,i) - s(mw,i)/6.d0 *
     &                     (1.d0 - (s(mw,i)*dtdxave)**2) * dq2

  110          continue
  120       continue
c
c
  140 continue
c
c     # update q by differencing correction fluxes 
c     ============================================
c
c     # (Note:  Godunov update has already been performed above)
c
      do m=1,meqn
	 do i=1,mx
	    q(m,i) = q(m,i) - dtdx(i) * (f(m,i+1) - f(m,i))
       enddo
      enddo
c
  900 continue
      return
      end

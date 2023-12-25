c =========================================================
      subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,
     &               wave,s,amdq,apdq)
c =========================================================
c
c     # solve Riemann problems for the 1D advection equation q_t + (u*q)_x = 0.
c
c       -----------------------------------------------------------
c     # In conservation form, with cell-centered velocities specified in
c     # the auxiliary variable
c     # aux(i,1)  =  u-velocity in cell i
c       -----------------------------------------------------------
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c     # On output, wave contains the waves,
c     #            s the speeds,
c     #            amdq the  left-going flux difference  A^- \Delta q
c     #            apdq the right-going flux difference  A^+ \Delta q
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routine step1, rp is called with ql = qr = q.
c
c
      implicit double precision (a-h,o-z)
      dimension   ql(meqn,1-mbc:maxmx+mbc)
      dimension   qr(meqn,1-mbc:maxmx+mbc)
      dimension auxl(maux,1-mbc:maxmx+mbc)
      dimension auxr(maux,1-mbc:maxmx+mbc)
      dimension    s(mwaves,1-mbc:maxmx+mbc)
      dimension wave(meqn,mwaves,1-mbc:maxmx+mbc)
      dimension amdq(meqn,1-mbc:maxmx+mbc)
      dimension apdq(meqn,1-mbc:maxmx+mbc)
c
      do i=2-mbc,mx+mbc
         ui  = auxl(1,i)
         uim = auxr(1,i-1)
         qi  = ql(1,i)
         qim = qr(1,i-1)
c
         if (uim .lt. 0.0d0 .and. ui .gt. 0.0d0) then
            wave(1,1,i) = qi - qim
            s(1,i) = 0.0d0
            amdq(1,i) = 0.0d0 - uim * qim
            apdq(1,i) = ui * qi - 0.0d0
         else if (uim .gt. 0.0d0 .and. ui .lt. 0.0d0) then
            stop "uim > 0 and ui < 0 not implemented"
         else if (ui .gt. 0.0d0) then
            qstar = uim*qim/ui
            wave(1,1,i) = qi - qstar
            s(1,i) = ui
            amdq(1,i) = 0.0d0
            apdq(1,i) = ui*qi - uim*qim
         else
            qstar = ui*qi/uim
            wave(1,1,i) = qstar - qim
            s(1,i) = uim
            amdq(1,i) = ui*qi - uim*qim
            apdq(1,i) = 0.0d0
         endif

      enddo
c
      return
      end

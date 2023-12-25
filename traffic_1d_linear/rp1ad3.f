c =========================================================
      subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,
     &               wave,s,amdq,apdq)
c =========================================================
c
c     # solve Riemann problems for the 1D advection equation q_t + (u*q)_x = 0.
c
c       -----------------------------------------------------------
c     # In conservation form, with interface velocities specified in
c     # the auxiliary variable
c     # aux(i,1)  =  u-velocity at left edge of cell i
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
c
c     # set fim1 = f(q_{i-1}) at left boundary for flux differencing below:
      i = 1-mbc
      fim1 = 0.5d0*(qr(1,i) + ql(1,i)) *
     &            (dmin1(auxl(1,i+1),0.0d0) + dmax1(auxl(1,i),0.0d0))
c
c
c
      do i=2-mbc,mx+mbc
c
c        # Compute the wave and speed
c
         u = auxl(1,i)
         wave(1,1,i) = ql(1,i) - qr(1,i-1)
         s(1,i) = u
c
c
c        # conservative form
c        -------------------
c        # amdq and apdq are chosen as flux differences for the
c        # conservative equation  q_t + (u*q)_x = 0
c
c        # compute the flux at the interface between cells i-1 and i:
         if (u.gt.0.d0) then
            f0 = u*qr(1,i-1)
         else
            f0 = u*ql(1,i)
         endif
c
c        # compute a value for the flux in cell i:
c        # note that we have velocities only at the interfaces
         fi = 0.5d0*(ql(1,i) + qr(1,i)) *
     &             (dmin1(auxl(1,i+1),0.0d0) + dmax1(u,0.0d0))
c
c        # flux differences:
         amdq(1,i) = f0 - fim1
         apdq(1,i) = fi - f0
c
         fim1 = fi

      enddo
c
      return
      end

c =========================================================
      subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,
     &               wave,s,amdq,apdq)
c =========================================================
c
c     # solve Riemann problems for the 1D Buckley-Leverett equation.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c     # On output, wave contains the waves, 
c     #            s the speeds, 
c     #            amdq the  left-going flux difference  A^- \Delta q
c     #            apdq the right-going flux difference  A^+ \Delta q
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routine step1, rp is called with ql = qr = q.
c
      implicit double precision (a-h,o-z)
      dimension   ql(meqn,1-mbc:maxmx+mbc)
      dimension   qr(meqn,1-mbc:maxmx+mbc)
      dimension    s(mwaves,1-mbc:maxmx+mbc)
      dimension wave(meqn,mwaves,1-mbc:maxmx+mbc)
      dimension amdq(meqn,1-mbc:maxmx+mbc)
      dimension apdq(meqn,1-mbc:maxmx+mbc)   
      common /comprob/ a
      
      do i=2-mbc,mx+mbc
         ur = ql(1,i)
         ul = qr(1,i-1)
         fr = ur**2/(ur**2 + a*(1.0d0-ur)**2)
         fl = ul**2/(ul**2 + a*(1.0d0-ul)**2)
c
         wave(1,1,i) = ur - ul
         if (ul.ne.ur) then
            s(1,i) = (fr-fl)/(ur-ul)
         else
            s(1,i) = 0.0d0
         endif

c        # compute left-going and right-going flux differences:
c        ------------------------------------------------------
c
         amdq(1,i) = 0.0d0
         apdq(1,i) = fr - fl
      enddo
c
      return
      end

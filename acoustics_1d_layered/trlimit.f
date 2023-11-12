c
c
c ===================================================================
      subroutine trlimit(meqn,mwaves,maux,mbc,mx,aux,wave,s,mthlim)
c ===================================================================
c
c     # Transmission-based limiter for acoustics equations
c     #    See Fogarty and LeVeque paper for a description.
c
c     # Use in place of limiter.f
c
c     # Assumes 
c     #    aux(i,1) = sound speed c in i'th cell
c     #    aux(i,2) = impedance Z in i'th cell
c     #    wave(i,*,1) is the left-going wave 
c     #      = alf_i^1 [Z_{i-1} ; 1]
c     #    wave(i,*,2) is the right-going wave 
c     #      = alf_i^2 [Z_i ; 1]
c     # 
c
c     --------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      dimension  aux(maux,1-mbc:mx+mbc)
      dimension    s(mwaves,1-mbc:mx+mbc)
      dimension wave(meqn,mwaves,1-mbc:mx+mbc)
      dimension mthlim(mwaves)
		external philim
c
c
c     # transmission-based limiter:
c     # use the fact that second component of wave is strength alpha since
c     # second component of eigenvector is 1.
c
       do 100 i=0,mx+1
c          # 1-wave at this cell and neighbor:
	   alf1i = wave(2,1,i)
	   if (alf1i .ne. 0.d0) then
	      alf1ip = wave(2,1,i+1)
c             # transmitted part of neighboring 1-wave:
	      alf1ipt = (2.d0*aux(1,i)/(aux(1,i-1)+aux(1,i))) * alf1ip
	      wlimitr = philim(alf1i, alf1ipt, mthlim(1))
	      do m=1,meqn
	         wave(m,1,i) = wlimitr * wave(m,1,i)
	         enddo
	      endif

c          # 2-wave at this cell and neighbor:
	   alf2i = wave(2,2,i)
	   if (alf2i .ne. 0.d0) then
	      alf2im = wave(2,2,i-1)
c             # transmitted part of neighboring 2-wave:
  	      alf2imt = (2.d0*aux(1,i-1)/(aux(1,i-1)+aux(1,i))) * alf2im
  	      wlimitr = philim(alf2i, alf2imt, mthlim(2))
	      do m=1,meqn
	         wave(m,2,i) = wlimitr * wave(m,2,i)
	         enddo
	      endif
  100      continue
c
      return
      end

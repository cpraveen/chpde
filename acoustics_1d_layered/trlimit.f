! ===================================================================
      subroutine trlimit (meqn, mwaves, maux, mbc, mx, aux, wave, s,
     &                    mthlim)
! ===================================================================
! 
!     # Transmission-based limiter for acoustics equations
!     #    See Fogarty and LeVeque paper for a description.
! 
!     # Use in place of limiter.f
! 
!     # Assumes
!     #    aux(i,1) = sound speed c in i'th cell
!     #    aux(i,2) = impedance Z in i'th cell
!     #    wave(i,*,1) is the left-going wave
!     #      = alf_i^1 [Z_{i-1} ; 1]
!     #    wave(i,*,2) is the right-going wave
!     #      = alf_i^2 [Z_i ; 1]
!     #
! 
!     ------------------------------------------------------------------
! 
      implicit double precision (a-h,o-z)
      dimension aux(maux,1-mbc:mx+mbc)
      dimension s(mwaves,1-mbc:mx+mbc)
      dimension wave(meqn,mwaves,1-mbc:mx+mbc)
      dimension mthlim(mwaves)
      external philim
! 
! 
!     # transmission-based limiter:
!     # use the fact that second component of wave is strength alpha sin
!     # second component of eigenvector is 1.
! 
      do i=0,mx+1
!        1-wave at this cell and neighbor:
         alf1i=wave(2,1,i)
         if (alf1i.ne.0.d0) then
            alf1ip=wave(2,1,i+1)
!           transmitted part of neighboring 1-wave:
            alf1ipt=(2.d0*aux(1,i)/(aux(1,i-1)+aux(1,i)))*alf1ip
            wlimitr=philim(alf1i,alf1ipt,mthlim(1))
            do m=1,meqn
               wave(m,1,i)=wlimitr*wave(m,1,i)
            end do
         end if
  
!        2-wave at this cell and neighbor:
         alf2i=wave(2,2,i)
         if (alf2i.ne.0.d0) then
            alf2im=wave(2,2,i-1)
!           transmitted part of neighboring 2-wave:
            alf2imt=(2.d0*aux(1,i-1)/(aux(1,i-1)+aux(1,i)))*alf2im
            wlimitr=philim(alf2i,alf2imt,mthlim(2))
            do m=1,meqn
               wave(m,2,i)=wlimitr*wave(m,2,i)
            end do
         end if
      end do
! 
      return
      end

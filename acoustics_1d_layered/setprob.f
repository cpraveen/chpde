      subroutine setprob
      implicit double precision (a-h,o-z)
      common /comaux/ rho(1000)
      common /comlim/ mylim
      common /combc/ omega
      common /comwall/ pi,t1,a1,tw1,t2,a2,tw2
  
!     # Set the material parameters for the acoustic equations
! 
      open (unit=7,file='setprob.data',status='old',form='formatted')
      open (unit=8,file='../setprob.rho',status='old',form='formatted')
  
      pi=4.d0*datan(1.d0)
  
!     Skip first 5 lines
      read (7,*)
      read (7,*)
      read (7,*)
      read (7,*)
      read (7,*)
  
!     # parameters for wall
      read (7,*) t1
      read (7,*) a1
      read (7,*) tw1
      read (7,*) t2
      read (7,*) a2
      read (7,*) tw2
  
!     # if mylim=1 then transmissve limiter is applied in rp1 rather tha
      read (7,*) mylim
  
!     # Piecewise constant medium with interfaces at x = 1, 2, ...
!     # Random densities can be created in matlab, e.g.:
!     #       rho = rand(1000,1) + 1;
!     #       save setprob.rho rho -ascii
  
!     # for periodic:
!     #       rho = mod(1:1000,2)*2 + 1;
  
      do i=1,1000
         read (8,*) rho(i)
      end do
! 
      return
      end

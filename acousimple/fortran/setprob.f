      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname

      common /cparam/ rho,bulk,cc,zz   
      common /cqinit/ beta
c
c     # Set the material parameters for the acoustic equations
c     # Passed to the Riemann solver rp1.f in a common block
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)


c     # density:
      read(iunit,*) rho

c     # bulk modulus:
      read(iunit,*) bulk

c     # sound speed:
      cc = dsqrt(bulk/rho)

c     # impedance:
      zz = cc*rho

c     # beta for initial conditions:
      read(iunit,*) beta

      close(iunit)

      return
      end




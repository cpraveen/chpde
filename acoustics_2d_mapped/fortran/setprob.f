      subroutine setprob
      implicit double precision (a-h,o-z)
      common /cparam/ rho,bulk,cc,zz

c
c     # Set the material parameters for the acoustic equations
c
      iunit = 7
      call opendatafile(iunit,'setprob.data')
c
c     # Density and sound speed:

      read(iunit,*) rho
      read(iunit,*) cc

      close(iunit)
c
c     # Compute bulk modulus and impedance:

      bulk = cc*cc*rho
      zz = rho*cc

      return
      end

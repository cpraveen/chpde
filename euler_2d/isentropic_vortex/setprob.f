c     ==================
      subroutine setprob
c     ==================

      implicit double precision (a-h,o-z)
      character*12 fname
      common /cparam/  gamma
      common /problem/ dmach, alpha, beta 
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                

c
       read(7,*) gamma
       read(7,*) dmach
       read(7,*) alpha
       read(7,*) beta

      return
      end

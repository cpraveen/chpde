      subroutine setprob
      implicit double precision (a-h,o-z)
c
c     # Copy this file to your directory and modify to set up problem
c     # parameters or read other data.
c
      common /sw/  g
      common/cdisc/ x0,y0,alf,beta,r0,idisc
      common /comic/ hl,hr,hul,hur

      g = 1.0d0

c     # data for flow into cylinder:
      idisc = 1
      x0 = -2.0d0
      y0 = 0.0d0
      alf = 1.0d0
      beta = 0.0d0

      hl = 4.0d0
      hr = 1.0d0
      hur = 0.0d0
      hul = hur + (hl-hr)*(hur/hr + dsqrt(g*hr + 0.5d0*g*(hl-hr)*
     &              (3.0d0 + (hl-hr)/hr)))

      return
      end

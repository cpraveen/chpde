import numpy as np
from numpy import pi,abs,cos,sin,tan,sqrt,fmin

r0 = 1.0
r1 = 5.0

def mapc2p(xc,yc):
    yc1 = fmin(yc, abs(pi/2-yc))
    yc1 = fmin(yc1, abs(pi-yc))
    yc1 = fmin(yc1, abs(3*pi/2-yc))
    yc1 = fmin(yc1, abs(2*pi-yc))
    r = r0 + (r1*sqrt(1.0 + tan(yc1)**2) - r0) * xc
    xp = r * cos(yc)
    yp = r * sin(yc)
    return xp, yp

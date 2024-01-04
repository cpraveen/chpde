import numpy as np

rin = 1.0         # radius of cylinder
rout = 50.0 * rin # radius of outer boundary

def mapc2p(xc,yc):
    r = rin * np.exp(xc * np.log(rout/rin))
    theta = 2.0 * np.pi * yc
    xp = r * np.cos(theta)
    yp = r * np.sin(theta)
    return xp, yp


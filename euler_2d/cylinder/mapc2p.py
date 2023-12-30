import numpy as np

rin = 1.0
rout = 50 * rin

def mapc2p(xc,yc):
    r = rin + xc * (rout - rin)
    theta = 2.0 * np.pi * yc
    xp = r * np.cos(theta)
    yp = r * np.sin(theta)
    return xp, yp


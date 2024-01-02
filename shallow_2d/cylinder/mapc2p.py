import numpy as np

def mapc2p(xc,yc):
    xp = xc * np.cos(yc)
    yp = xc * np.sin(yc)
    return xp, yp

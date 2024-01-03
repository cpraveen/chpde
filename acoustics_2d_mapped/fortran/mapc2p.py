from numpy import abs

def mapc2p(xc,yc):
    xp = xc + (abs(yc+0.2) + 0.8)/2.0
    yp = yc
    return xp, yp

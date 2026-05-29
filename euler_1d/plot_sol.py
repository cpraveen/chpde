import numpy as np
import matplotlib.pyplot as plt
import clawpack.pyclaw as pyclaw
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-frame', type=int, help='Frame number', default=10)
parser.add_argument('-gamma', type=float, help='Gas gamma', default=1.4)
args = parser.parse_args()

gamma = args.gamma
frame = args.frame

f  = pyclaw.Solution()
f.read(frame)

t  = f.state.t
print("Solution at time t = ", t)

x   = f.state.grid.x.centers
rho = f.state.q[0,:]
v   = f.state.q[1,:] / rho
E   = f.state.q[2,:]
p   = (gamma - 1.0) * (E - 0.5 * rho * v**2)

plt.figure()
plt.plot(x, rho, '-')
plt.ylabel("Density")

plt.figure()
plt.plot(x, v, '-')
plt.ylabel("Velocity")

plt.figure()
plt.plot(x, p, '-')
plt.ylabel("Pressure")

plt.show()

#!/usr/bin/env python
# encoding: utf-8

"""
2D Euler Riemann problem
================================

Solve the Euler equations of compressible fluid dynamics:

.. math::
    \rho_t + (\rho u)_x + (\rho v)_y & = 0 \\
    (\rho u)_t + (\rho u^2 + p)_x + (\rho uv)_y & = 0 \\
    (\rho v)_t + (\rho uv)_x + (\rho v^2 + p)_y & = 0 \\
    E_t + (u (E + p) )_x + (v (E + p))_y & = 0.

Here :math:`\rho` is the density, (u,v) is the velocity, and E is the total energy.
The initial condition is one of the 2D Riemann problems from the paper of
Liska and Wendroff, Case 6 in Table 4.3, see
    https://doi.org/10.1137/S1064827502402120

This example shows how to use VisClaw's Iplot class for simple interactive plotting.
"""
from clawpack import pyclaw
from clawpack import riemann
from clawpack.riemann.euler_4wave_2D_constants import density, x_momentum, \
        y_momentum, energy, num_eqn

def load_frame(frame_number):
    from clawpack.pyclaw import Solution

    try:
        return Solution(frame_number)
    except:
        print("Reached end of frames")
        exit()

def plot_frame(frame):
    import matplotlib.pyplot as plt
    q = frame.q
    x, y = frame.state.grid.c_centers
    plt.cla()
    #plt.pcolormesh(x, y, q[density,...], cmap='viridis')
    plt.contour(x, y, q[density,...], levels=20)
    plt.title('Density, t=' + str(frame.t))
    plt.axis('equal')

def plot_results():
    from clawpack.visclaw import iplot
    ip = iplot.Iplot(load_frame,plot_frame)
    ip.plotloop()

solver = pyclaw.ClawSolver2D(riemann.euler_4wave_2D)
solver.all_bcs = pyclaw.BC.extrap

mx, my = 200, 200
domain = pyclaw.Domain([0.0,0.0],[1.0,1.0],[mx,my])
solution = pyclaw.Solution(num_eqn,domain)
gamma = 1.4
solution.problem_data['gamma']  = gamma

# Set initial data
xx,yy = domain.grid.p_centers
l = xx<0.5; r = xx>=0.5; b = yy<0.5; t = yy>=0.5
u = 0.75*t - 0.75*b
v = 0.5*l  - 0.5*r
solution.q[density,...] = 2.0*l*t + 1.0*l*b + 1.0*r*t + 3.0*r*b
solution.q[x_momentum,...] = solution.q[density,...] * u
solution.q[y_momentum,...] = solution.q[density,...] * v
solution.q[energy,...] = 0.5*solution.q[density,...]*(u**2 + v**2) \
                         + 1.0/(gamma-1.0)

claw = pyclaw.Controller()
claw.tfinal = 0.3
claw.solution = solution
claw.solver = solver

status = claw.run()

plot_results()

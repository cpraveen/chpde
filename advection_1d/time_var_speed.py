#!/usr/bin/env python
# encoding: utf-8

r"""
One-dimensional advection
=========================

Solve the linear advection equation:

.. math:: 
    q_t + u(t) q_x = 0.

Here q is the density of some conserved quantity and u is the velocity.

The initial condition is a Gaussian and the boundary conditions are periodic.
"""
from __future__ import absolute_import
import numpy as np
from clawpack import riemann

omega = 2.0 * np.pi
def velocity(t):
    return np.cos(omega * t)

def exact(x,t):
    beta = 100; gamma = 0; x0 = 0.5
    xc = x - np.sin(omega * t) / omega
    return np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))

def plot_exact(current_data):
    from pylab import plot,legend
    t = current_data.t
    x = np.linspace(0,1,100)
    plot(x, exact(x,t),'r-',label='Exact')
    legend()

def b4step(solver,state):
    t = state.t
    dt = solver.dt
    state.problem_data['u'] = velocity(t + 0.5*dt)

def setup(nx=100, kernel_language='Python', solver_type='classic',
          outdir='./_output'):

    from clawpack import pyclaw

    if kernel_language == 'Fortran':
        riemann_solver = riemann.advection_1D
    elif kernel_language == 'Python':
        riemann_solver = riemann.advection_1D_py.advection_1D
            
    solver = pyclaw.ClawSolver1D(riemann_solver)

    solver.kernel_language = kernel_language
    solver.limiters = 0
    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.before_step = b4step
    solver.cfl_max = 1.0
    solver.cfl_desired = 0.9

    x = pyclaw.Dimension(0.0,1.0,nx,name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,solver.num_eqn)

    state.problem_data['u'] = 0.0  # Advection velocity

    # Initial data
    xc = state.grid.x.centers
    state.q[0,:] = exact(xc, 0.0)

    claw = pyclaw.Controller()
    claw.keep_copy = False
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.tfinal = 1.0
    claw.output_style = 1
    claw.num_output_times = 20
    claw.setplot = setplot

    return claw

def setplot(plotdata):
    """ 
    Plot solution using VisClaw.
    """ 
    plotdata.clearfigures()  # clear any old figures,axes,items data

    plotfigure = plotdata.new_plotfigure(name='q', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.ylimits = [-0.2,1.0]
    plotaxes.title = 'q'
    plotaxes.afteraxes = plot_exact

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'o'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':2,'markersize':5}
    
    return plotdata

 
if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)

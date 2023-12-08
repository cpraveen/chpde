#!/usr/bin/env python
# encoding: utf-8

r"""
Shallow water flow
==================

Solve the one-dimensional shallow water equations with trace:

.. math::
    h_t + (hu)_x & = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x & = 0 \\
    c + u c_x = 0

Here h is the depth, u is the velocity, c is tracer concentration, and g is the 
gravitational constant.  The default initial condition used here models a 
dam break.
"""

from __future__ import absolute_import
import numpy as np
from clawpack import riemann
from clawpack.riemann.shallow_roe_tracer_1D_constants import depth, momentum, tracer, num_eqn

def setup(IC='dam-break',use_petsc=False,
          outdir='./_output',solver_type='classic',
          disable_output=False):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    riemann_solver = riemann.shallow_roe_tracer_1D
 
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver1D(riemann_solver)
        solver.limiters = pyclaw.limiters.tvd.vanleer
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D(riemann_solver)

    solver.kernel_language = 'Fortran'

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap

    xlower, xupper = -5.0, 5.0
    mx = 500
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,num_eqn)

    # Gravitational constant
    state.problem_data['grav'] = 1.0
    state.problem_data['dry_tolerance'] = 1e-3
    state.problem_data['sea_level'] = 0.0
    
    xc = state.grid.x.centers


    if IC=='dam-break':
        tf = 2.0
        hl = 3.0
        ul = 0.0
        hr = 1.0
        ur = 0.0
        x0 = 0.0
        state.q[depth,:] = hl * (xc <= x0) + hr * (xc > x0)
        state.q[momentum,:] = hl*ul * (xc <= x0) + hr*ur * (xc > x0)
    elif IC=='2-shock':
        tf = 2.0
        hl = 1.0
        ul = 1.0
        hr = 1.0
        ur = -1.0
        x0 = 0.0
        state.q[depth,:] = hl * (xc <= x0) + hr * (xc > x0)
        state.q[momentum,:] = hl*ul * (xc <= x0) + hr*ur * (xc > x0)
    elif IC=='2-rare':
        tf = 3.0
        hl = 1.0
        ul = -0.5
        hr = 1.0
        ur = 0.5
        x0 = 0.0
        state.q[depth,:] = hl * (xc <= x0) + hr * (xc > x0)
        state.q[momentum,:] = hl*ul * (xc <= x0) + hr*ur * (xc > x0)
    elif IC=='perturb':
        tf = 3.0
        eps = 0.1
        x0 = 0.0
        state.q[depth,:] = 1.0 + eps*np.exp(-(xc-x0)**2/0.5)
        state.q[momentum,:] = 0.0

    state.q[tracer,:] = (xc <= x0)*(1.0) + (xc > x0) * (0.0)

    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.tfinal = tf
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot

    return claw


#--------------------------
def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for depth
    plotfigure = plotdata.new_plotfigure(name='Water height', figno=0)
    plotfigure.kwargs = {'figsize':[8,10]}

    # Figure for Water height
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-5.0,5.0]
    plotaxes.title = 'Water height'
    plotaxes.axescmd = 'subplot(311)'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = depth
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}

    # Figure for momentum[1]
    #plotfigure = plotdata.new_plotfigure(name='Momentum', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(312)'
    plotaxes.xlimits = [-5.0,5.0]
    plotaxes.title = 'Momentum'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = momentum
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}
    
    # Figure for tracer
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = [-5.0,5.0]
    plotaxes.title = 'Tracer'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = tracer
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}

    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)

#!/usr/bin/env python
# encoding: utf-8

r"""
Shallow water flow
==================

Solve the one-dimensional shallow water equations including bathymetry:

.. math::
    h_t + (hu)_x & = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x & = -g h b_x.

Here h is the depth, u is the velocity, g is the gravitational constant, and b
the bathymetry.  
"""

from __future__ import absolute_import
import numpy as np
#from clawpack import riemann
import shallow_1d

vmax = 5.0
x0 = 300.0 * np.sqrt(10.0)
H0 = 10.0
g = 9.81
omega = np.sqrt(2.0 * g)/x0
dry_tolerance = 1.0e-2
xlower, xupper = -5000.0, 5000.0

def eta(x,t):
    return (H0 - vmax**2/(4.0*g) * (1.0 + np.cos(2*omega*t)) 
            - np.sqrt(2/g) * vmax/x0 * np.cos(omega*t) * x)

def setup(kernel_language='Fortran', solver_type='classic', use_petsc=False,
          outdir='./_output'):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    solver = pyclaw.ClawSolver1D(shallow_1d)
    solver.kernel_language = 'Fortran'
    solver.order = 2
    solver.limiters = pyclaw.limiters.tvd.MC
    solver.fwave = True
    solver.num_waves = 2
    solver.num_eqn = 2
    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap

    x = pyclaw.Dimension(xlower, xupper, 500, name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain, 2, 1)

    # Gravitational constant
    state.problem_data['grav'] = g
    state.problem_data['dry_tolerance'] = dry_tolerance
    state.problem_data['sea_level'] = 0.0

    xc = state.grid.x.centers
    state.aux[0, :] = (xc/x0)**2
    state.q[0, :] = np.maximum(0.0, eta(xc, 0.0) - state.aux[0,:])
    state.q[1, :] = 0.0

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.tfinal = 5.0e3
    claw.cfl_max = 1.0
    claw.cfl_desired = 0.8
    claw.num_output_times = 20
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.setplot = setplot
    claw.write_aux_init = True

    if outdir is not None:
        claw.outdir = outdir
    else:
        claw.output_format = None


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

    # Plot variables
    def bathy(current_data):
        return current_data.aux[0, :]

    def eta(current_data):
        return current_data.q[0, :] + bathy(current_data)

    def velocity(current_data):
        h = current_data.q[0, :]
        m = current_data.q[1, :]
        v = ((h < dry_tolerance) * 0.0 + 
            (h > dry_tolerance) * m / h)
        return v

    rgb_converter = lambda triple: [float(rgb) / 255.0 for rgb in triple]

    # Figure for depth
    plotfigure = plotdata.new_plotfigure(name='Depth', figno=0)
    plotfigure.kwargs = {'layout': 'tight'}

    # Axes for water depth
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [xlower, xupper]
    #plotaxes.ylimits = [-1.1, 0.2]
    plotaxes.title = 'Water Depth'
    plotaxes.axescmd = 'subplot(211)'

    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = eta
    plotitem.plot_var2 = bathy
    plotitem.color = rgb_converter((67,183,219))

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = bathy
    plotitem.color = 'k'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = eta
    plotitem.color = 'k'

    # Axes for velocity
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = [xlower, xupper]
    #plotaxes.ylimits = [-0.5, 0.5]
    plotaxes.title = 'Velocity'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = velocity
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}
    
    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)

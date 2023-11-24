#!/usr/bin/env python
# encoding: utf-8
    
from __future__ import absolute_import
def setup(iplot=False,htmlplot=True,outdir='./_output'):
    r"""Produces the output shown in Figures 3.1 and 3.8 of the FVM book.
    These involve simple waves in the acoustics system."""
    from clawpack import pyclaw
    from clawpack import riemann
    import numpy as np

    solver = pyclaw.ClawSolver1D(riemann.acoustics_1D)

    solver.limiters = pyclaw.limiters.tvd.MC
    solver.bc_lower[0] = pyclaw.BC.wall
    solver.bc_upper[0] = pyclaw.BC.extrap

    x = pyclaw.Dimension(-1.0,1.0,800,name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain, solver.num_eqn)

    # Set problem-specific variables
    rho = 1.0
    bulk = 0.25
    state.problem_data['rho']=rho
    state.problem_data['bulk']=bulk
    state.problem_data['zz']=np.sqrt(rho*bulk)
    state.problem_data['cc']=np.sqrt(bulk/rho)

    # Set the initial condition
    xc = domain.grid.x.centers
    state.q[0,:] = 0.5*np.exp(-80 * xc**2) + 0.5*(np.abs(xc+0.2)<0.1)
    state.q[1,:] = 0.0
    
    # Set up the controller
    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.iplot = iplot
    claw.htmlplot = htmlplot
    claw.outdir = outdir
    claw.tfinal = 3.0
    claw.output_format = 'ascii'
    claw.num_output_times   = 30
    claw.setplot = setplot

    return claw


""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
"""
# --------------------------
def setplot(plotdata):
    # --------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Acoustics', figno=1)
    plotfigure.kwargs = {'figsize' : [12, 12]}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-1.0, 1.0]
    plotaxes.title = 'Pressure'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    plotitem.show = True       # show on plot?
    plotitem.kwargs = {'linewidth': 2, 'markersize': 4}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-1.0, 1.0]
    plotaxes.title = 'Velocity'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 1
    plotitem.plotstyle = '-o'
    plotitem.color = 'r'
    plotitem.show = True       # show on plot?
    plotitem.kwargs = {'linewidth': 2, 'markersize': 4}

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup, setplot)

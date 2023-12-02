#!/usr/bin/env python
# encoding: utf-8
r"""
Advection with div-free rotating velocity
==============================

Solve the linear non-conservative advection equation:

.. math::
    q_t + u(x,y) q_x + v(x,y) q_y = 0

in square domain.

Here q is the density of some quantity and (u,v) is the velocity
field.  We take a rotational velocity field: :math:`u = 2y, v = -2x`.
"""
import numpy as np

def qinit(state):
    x, y = state.grid.p_centers
    r = np.sqrt((x + 0.45)**2 + y**2)
    state.q[0,:,:]  = (x > 0.1)*(x < 0.6)*(y > -0.25)*(y < 0.25)*(1.0)
    state.q[0,:,:] += (r < 0.35)*(1.0 - r/0.35)

def psi(Xp,Yp):
    """ 
    Calculates the stream function in physical space.
    Clockwise rotation. One full rotation corresponds to pi (second).
    """
    return Xp**2 + Yp**2


def setup(use_petsc=False,outdir='./_output',solver_type='classic'):
    from clawpack import riemann
    from clawpack.riemann.vc_advection_2D_constants import num_aux

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    # Riemann solver
    rsolver = riemann.vc_advection_2D

    # Main solver
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D(rsolver)
        solver.dimensional_split = False
        solver.transverse_waves = 2
        solver.order = 2
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(rsolver)

    solver.all_bcs = pyclaw.BC.extrap
    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[1] = pyclaw.BC.extrap
    solver.aux_bc_upper[1] = pyclaw.BC.extrap
    solver.cfl_max = 1.0
    solver.cfl_desired = 0.9
    solver.limiters = pyclaw.limiters.tvd.MC

    xlower, xupper = -1.0, 1.0
    ylower, yupper = -1.0, 1.0
    mx, my = 80, 80

    x     = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y     = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = pyclaw.Domain([x,y])
    domain.grid.num_ghost = solver.num_ghost

    state = pyclaw.State(domain,solver.num_eqn,num_aux)

    # Initial condition
    qinit(state)

    # Edge velocities
    dx, dy = state.grid.delta
    Xe, Ye = domain.grid.p_nodes
    # u(x_(i-1/2),y_j)
    state.aux[0,:,:] =  ( psi(Xe[:-1,1:],Ye[:-1,1:]) - psi(Xe[:-1,:-1],Ye[:-1,:-1]) ) / dy
    # v(x_i,y_(j-1/2))
    state.aux[1,:,:] = -( psi(Xe[1:,:-1],Ye[1:,:-1]) - psi(Xe[:-1,:-1],Ye[:-1,:-1]) ) / dx

    claw = pyclaw.Controller()
    claw.tfinal = np.pi
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot
    claw.keep_copy = True

    return claw

def setplot(plotdata):
    """ 
    Plot solution using VisClaw.
    """ 
    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for pcolor plot
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'q[0]'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 1.0
    plotitem.add_colorbar = True
    
    # Figure for contour plot
    plotfigure = plotdata.new_plotfigure(name='contour', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'q[0]'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = 0
    plotitem.contour_nlevels = 20
    plotitem.contour_min = 0.01
    plotitem.contour_max = 0.99
    plotitem.amr_contour_colors = ['b','k','r']
    
    return plotdata

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)

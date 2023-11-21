#!/usr/bin/env python
# encoding: utf-8
r"""
Woodward-Colella blast wave problem
===================================

Solve the one-dimensional Euler equations for inviscid, compressible flow:

.. math::
    \rho_t + (\rho u)_x & = 0 \\
    (\rho u)_t + (\rho u^2 + p)_x & = 0 \\
    E_t + (u (E + p) )_x & = 0.

The fluid is an ideal gas, with pressure given by :math:`p=\rho (\gamma-1)e` where
e is internal energy.

This script runs the Woodward-Colella blast wave interaction problem,
involving the collision of two shock waves.

This example also demonstrates:

 - How to use a total fluctuation solver in SharpClaw
 - How to use characteristic decomposition with an evec() routine in SharpClaw
"""
from clawpack import riemann
from clawpack.riemann.euler_with_efix_1D_constants import *

gamma = 1.4 # Ratio of specific heats

def setup(order=2,flux='roe',use_petsc=False,outdir='./_output',
          solver_type='classic',kernel_language='Fortran',tfluct_solver=False):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if kernel_language =='Python':
        if flux == 'hll' or flux == 'hlle':
            rs = riemann.euler_1D_py.euler_hll_1D
        elif flux == 'hllc':
            rs = riemann.euler_1D_py.euler_hllc_1D
        elif flux == 'roe':
            rs = riemann.euler_1D_py.euler_roe_1D
        else:
            print("Python kernel:")
            print("   Use flux = hll, hllc, roe")
            exit(1)
    elif kernel_language =='Fortran':
        if flux == 'hll' or flux == 'hlle':
            rs = riemann.euler_hlle_1D
        elif flux == 'roe':
            # This always has efix enabled
            rs = riemann.euler_with_efix_1D
        else:
            print("Fortran kernel:")
            print("   Use flux = hll, roe")
            exit(1)

    if solver_type=='sharpclaw':
        assert flux == 'roe'
        solver = pyclaw.SharpClawSolver1D(rs)
        solver.time_integrator = 'SSP33'
        solver.cfl_max = 0.65
        solver.cfl_desired = 0.6
        solver.tfluct_solver = tfluct_solver
        if solver.tfluct_solver:
            import euler_tfluct
            solver.tfluct = euler_tfluct
        solver.lim_type = 1 # 0=none, 1=TVD, 2=weno
        solver.limiters = pyclaw.limiters.tvd.minmod
        solver.char_decomp = 2
        import euler_sharpclaw1
        solver.fmod = euler_sharpclaw1
    elif solver_type=='classic':
        solver = pyclaw.ClawSolver1D(rs)
        solver.order = order
        solver.limiters = pyclaw.limiters.tvd.MC

    solver.kernel_language = kernel_language

    solver.bc_lower[0]=pyclaw.BC.wall
    solver.bc_upper[0]=pyclaw.BC.wall

    mx = 800
    x = pyclaw.Dimension(0.0,1.0,mx,name='x')
    domain = pyclaw.Domain([x])
    state = pyclaw.State(domain,num_eqn)

    state.problem_data['gamma'] = gamma
    if kernel_language =='Python':
        state.problem_data['efix'] = False

    x = state.grid.x.centers

    pre = (x < 0.1)*1.0e3 + (0.1 <= x)*(x < 0.9)*1.0e-2 + (0.9 <= x)*1.0e2
    state.q[density ,:] = 1.0
    state.q[momentum,:] = 0.0
    state.q[energy  ,:] = pre / (gamma - 1.0)

    claw = pyclaw.Controller()
    claw.tfinal = 0.038
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.num_output_times = 10
    claw.outdir = outdir
    claw.setplot = setplot
    claw.keep_copy = True

    return claw

#--------------------------
def velocity(current_data):
    q  = current_data.q
    return q[1,:]/q[0,:]

#--------------------------
def pressure(current_data):
    q  = current_data.q
    rho = q[0,:]
    vel = q[1,:] / rho
    return  (gamma - 1.0) * (q[2,:] - 0.5 * rho * vel**2)

#--------------------------
def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    plotdata.clearfigures()  # clear any old figures,axes,items data

    plotfigure = plotdata.new_plotfigure(name='', figno=0)
    plotfigure.kwargs = {'figsize': [8,10], 'layout': 'tight'}

    # Density
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(311)'
    plotaxes.title = 'Density'

    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = density
    plotitem.kwargs = {'linewidth':2}
    
    # Velocity
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(312)'
    plotaxes.title = 'Velocity'

    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = velocity
    plotitem.kwargs = {'linewidth':2}
    
    # Pressure
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.title = 'Pressure'

    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = pressure
    plotitem.kwargs = {'linewidth':2}

    return plotdata

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)

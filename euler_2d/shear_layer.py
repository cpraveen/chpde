'''
G. Leidi, R. Andrassy, W. Barsukow, J. Higl, P. V. F. Edelmann, F. K. Röpke
Performance of high-order Godunov-type methods in simulations of astrophysical low Mach number flows
https://arxiv.org/abs/2402.16706
See Section 3.1: Kelvin-Helmholtz instability
'''

from clawpack import riemann
from clawpack.riemann.euler_4wave_2D_constants import density, x_momentum, \
        y_momentum, energy, num_eqn
from clawpack.visclaw import colormaps
import numpy as np

def setplot(plotdata):
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for density - pcolor
    plotfigure = plotdata.new_plotfigure(name='Density', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.scaled = True
    plotaxes.title = 'Density'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = density
    plotitem.contour_nlevels = 20
    #plotitem.contour_min = 0.5
    #plotitem.contour_max = 1.0

    return plotdata


def eta(y):
    X1 = 0.5 * (1.0 + np.sin(16 * np.pi * (y + 0.25)));
    X2 = 0.5 * (1.0 - np.sin(16 * np.pi * (y - 0.25)))
    return (y > -9/32)*(y < -7/32)*X1                  \
            + (y >= -7/32)*(y <= 7/32)*np.ones_like(y) \
            + (y > 7/32)*(y < 9/32)*X2                 \
            + 0.0

def setup(use_petsc=False,riemann_solver='roe'):
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if riemann_solver.lower() == 'roe':
        solver = pyclaw.ClawSolver2D(riemann.euler_4wave_2D)
        solver.transverse_waves = 2
    elif riemann_solver.lower() == 'hlle':
        solver = pyclaw.ClawSolver2D(riemann.euler_hlle_2D)
        solver.transverse_waves = 0
        solver.cfl_desired = 0.4
        solver.cfl_max = 0.5
    solver.limiters = pyclaw.limiters.tvd.MC
    solver.all_bcs = pyclaw.BC.periodic

    mx, my = 2*128, 2*64
    domain = pyclaw.Domain([0.0,-0.5],[2.0,0.5],[mx,my])
    solution = pyclaw.Solution(num_eqn,domain)
    gamma = 1.4
    solution.problem_data['gamma']  = gamma

    # Set initial data
    # Initial level of mach number
    M0 = 0.1

    x, y = domain.grid.p_centers
    rho  = gamma + 1.0e-3 * (1.0 - 2.0 * eta(y))
    u    = M0 * (1.0 - 2.0 * eta(y))
    v    = 0.1 * M0 * np.sin(2.0 * np.pi * x)
    p    = np.ones_like(rho)

    solution.q[density,...] = rho
    solution.q[x_momentum,...] = rho * u
    solution.q[y_momentum,...] = rho * v
    solution.q[energy,...] = 0.5 * rho * (u**2 + v**2) + p / (gamma - 1.0)

    claw = pyclaw.Controller()
    claw.tfinal = 0.8 / M0
    claw.num_output_times = 80
    claw.solution = solution
    claw.solver = solver

    claw.output_format = 'ascii'    
    claw.outdir = "./_output"
    claw.setplot = setplot

    return claw

if __name__ == "__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup, setplot)

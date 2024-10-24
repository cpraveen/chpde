from clawpack import riemann
from clawpack.riemann.euler_4wave_2D_constants import density, x_momentum, \
        y_momentum, energy, num_eqn
from clawpack.visclaw import colormaps
import numpy as np

def setplot(plotdata):
    r"""Plotting settings
    Should plot two figures both of density.
    """

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
    plotitem.contour_min = 0.5
    plotitem.contour_max = 1.0

    # Figure for density - Schlieren
    plotfigure = plotdata.new_plotfigure(name='Schlieren', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Schlieren'
    plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_schlieren')
    plotitem.schlieren_cmin = 0.0
    plotitem.schlieren_cmax = 1.0
    plotitem.plot_var = density
    plotitem.add_colorbar = False
    
    return plotdata


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

    L, mx, my = 5.0, 100, 100
    domain = pyclaw.Domain([-L,-L],[L,L],[mx,my])
    solution = pyclaw.Solution(num_eqn,domain)
    gamma = 1.4
    solution.problem_data['gamma']  = gamma

    # Set initial data
    mach, alpha, beta = 0.5, 45.0, 5.0
    vx0 = mach * np.cos(alpha*np.pi/180)
    vy0 = mach * np.sin(alpha*np.pi/180)
    period = np.sqrt(2.0)*(2.0*L) / mach

    x, y = domain.grid.p_centers
    r2 = x**2 + y**2
    T = 1.0 - (gamma-1.0) * beta**2 / (8.0 * gamma * np.pi**2) * np.exp(1.0-r2)
    rho = T**(1.0/(gamma-1.0))
    u = vx0 - beta/(2.0*np.pi)*y*np.exp(0.5*(1.0-r2))
    v = vy0 + beta/(2.0*np.pi)*x*np.exp(0.5*(1.0-r2))
    p = rho**gamma

    solution.q[density,...] = rho
    solution.q[x_momentum,...] = rho * u
    solution.q[y_momentum,...] = rho * v
    solution.q[energy,...] = 0.5 * rho * (u**2 + v**2) + p / (gamma - 1.0)

    claw = pyclaw.Controller()
    claw.tfinal = period
    claw.num_output_times = 40
    claw.solution = solution
    claw.solver = solver

    claw.output_format = 'ascii'    
    claw.outdir = "./_output"
    claw.setplot = setplot

    return claw

if __name__ == "__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup, setplot)

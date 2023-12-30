
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 
from __future__ import absolute_import
import numpy as np

#--------------------------
def setplot(plotdata=None):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    from clawpack.visclaw import colormaps

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Compute x velocity
    def v1(current_data):
        q  = current_data.q
        h  = q[0,:,:]
        v1 = q[1,:,:]/h
        return v1

    # Figure for pcolor of average depth
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Pcolor_h', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-5,5]
    plotaxes.ylimits = [-5,5]
    plotaxes.title = 'h'
    plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.blue_yellow_red
    #plotitem.pcolor_cmin = 0.0
    #plotitem.pcolor_cmax = 0.02
    plotitem.add_colorbar = True

    # Figure for pcolor of v1 
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Pcolor_v1', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-5,5]
    plotaxes.ylimits = [-5,5]
    plotaxes.title = 'v1'
    plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = v1
    plotitem.pcolor_cmap = colormaps.blue_yellow_red
    #plotitem.pcolor_cmin = 0.0
    #plotitem.pcolor_cmax = 1.0
    plotitem.add_colorbar = True

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.html_movie = 'JSAnimation'      # new style, or "4.x" for old style
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

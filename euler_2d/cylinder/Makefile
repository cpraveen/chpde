
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = classic                  # Clawpack package to use
EXE = xclaw                         # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?=  -O3 -Wall

# ---------------------------------
# package sources for this program:
# ---------------------------------

include $(CLAW)/classic/src/2d/Makefile.classic_2d

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \


# ---------------------------------
# List of custom sources for this program:
# ---------------------------------

MODULES = \

SOURCES = \
  $(CLAW)/riemann/src/rpn2_euler_mapgrid.f90 \
  $(CLAW)/riemann/src/rpt2_euler_mapgrid.f90 \
  $(CLAW)/riemann/src/euler_roe_solver_mapgrid.f90 \
  $(CLAW)/riemann/src/getquadinfo_mapgrid.f90 \
  qinit.f90 \
  setprob.f90 \
  setaux.f90 \
  mapc2p.f90 \
  bc2.f

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import clawpack.pyclaw as pyclaw
import argparse
from acoustics_2d_inclusions import inclusion_mapping

from matplotlib import rcParams
rcParams['font.size'] = 12
rcParams['font.family'] = 'serif'
rcParams['figure.autolayout'] = True
rcParams['lines.linewidth'] = 1
rcParams['lines.markersize'] = 6
rcParams['axes.titlesize'] = 12
rcParams['axes.labelsize'] = 12

#------------------------------------------------------------------------------
def read(frame, dir, petsc):
    f = pyclaw.Solution()
    if petsc:
        f.read(frame, dir, read_aux=False, file_prefix='claw', 
               file_format='petsc')
    else:
        f.read(frame, dir, read_aux=False)
    return f
#------------------------------------------------------------------------------
# Data Reading

# Defaults
dir = "_output"

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-d', help='Directory', default=dir)
parser.add_argument('-petsc', help='Read PETSc output', action='store_true')
args = parser.parse_args()
dir = args.d

print('Dir   = ', dir)

# Read first frame and make grid
f = read(0, dir, args.petsc)
f.state.grid.mapc2p = inclusion_mapping
x,y  = f.state.grid.p_nodes
nx, ny = f.state.grid.num_cells

plt.figure(layout='tight')
for i in range(nx+1):
    plt.plot(x[i,:],y[i,:],'k-')
# Plot radial lines
for j in range(ny+1):
    plt.plot(x[:,j],y[:,j],'k-')
plt.title('Contours of q \n')
plt.xlabel('x'); plt.ylabel('y')
plt.axis('equal')
plt.show()

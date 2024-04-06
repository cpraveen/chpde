from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import clawpack.pyclaw as pyclaw
import argparse

from matplotlib import rcParams
rcParams['font.size'] = 12
rcParams['font.family'] = 'serif'
rcParams['figure.autolayout'] = True
rcParams['lines.linewidth'] = 2
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
def box():
    x = [-1,1,1,-1,-1]
    y = [-1, -1, 1, 1, -1]
    plt.plot(x,y,'-k')
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
x, y  = f.state.grid.p_centers
nx, ny = f.state.grid.num_cells

# Contour levels to draw
levels = np.arange(0.05, 1.0, 0.1)

plt.figure(layout='tight')
plt.title('Contours of q \n')
plt.xlabel('x'); plt.ylabel('y')
plt.box(False); plt.axis('off')
plt.show(block=False)

for frame in range(1000):
    try:
        f = read(frame, dir, args.petsc)
    except:
        print("Reached end of frames")
        input("Press Enter to exit")
        break

    t  = f.state.t

    # 2d plot data
    q  = f.state.q[0,:,:]
    #---------------------------------------------------------------------------
    print('Frame, t = ', frame, t)
    t = "{:.3f}".format(t)

    #plt.contourf(x,y,q,levels=50,cmap='jet')
    C = plt.contour(x,y,q,levels=levels,linewidths=1)
    plt.title('Time ='+str(t), loc='right')
    box(); plt.axis('equal'); plt.draw()
    input("Press Enter to continue...")
    C.remove()

#------------------------------------------------------------------------------

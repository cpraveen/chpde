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
levels = np.arange(0.0, 2.0, 0.1)

fig = plt.figure(layout='tight')
ax = fig.add_subplot(111, projection='3d')
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
    h  = f.state.q[0,:,:]

    #---------------------------------------------------------------------------
    print('Frame, t = ', frame, t)
    t = "{:.3f}".format(t)

    surf = ax.plot_surface(x,y,h,color='silver',shade=True,linewidths=1)
    ax.set_title('Time ='+str(t), loc='right')
    ax.set_zlabel('h')
    ax.set_zlim(0.0, 2.0)
    cont = ax.contour(x,y,h,offset=0.01,zdir='z',levels=20)
    plt.draw()
    input("Press Enter to continue...")
    surf.remove(); cont.remove()

#------------------------------------------------------------------------------

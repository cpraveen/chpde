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
def read(frame):
    f = pyclaw.Solution()
    if args.petsc:
        f.read(frame, dir, read_aux=False, file_prefix='claw', 
               file_format='petsc')
    else:
        f.read(frame, dir, read_aux=False)
    return f
#------------------------------------------------------------------------------
# Data Reading

# Defaults
frame = 10
dir = "_output"

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', type=int, help='Frame number', default=frame)
parser.add_argument('-d', help='Directory', default=dir)
parser.add_argument('-petsc', help='Read PETSc output', action='store_true')
args = parser.parse_args()
dir = args.d

print('Dir   = ', dir)

# Read first frame and make grid
f = read(0)
x  = f.state.grid.x.centers
y  = f.state.grid.y.centers
nx, ny = len(x), len(y)
x, y = np.meshgrid(x, y, indexing="ij")

# Contour levels to draw
levels = np.arange(0.05, 1.0, 0.1)

plt.figure(layout='tight')
plt.title('Contours of q \n')
plt.xlabel('x')
plt.ylabel('y')
plt.show(block=False)

for frame in range(1000):
    try:
        f = read(frame)
    except:
        print("End")
        break

    t  = f.state.t

    # 2d plot data
    q  = f.state.q[0,:,:]
    #------------------------------------------------------------------------------
    print('Frame, t = ', frame, t)
    t = "{:.3f}".format(t)

    #plt.contourf(x,y,q,levels=50,cmap='jet')
    C = plt.contour(x,y,q,levels=levels,linewidths=1)
    plt.title('Time ='+str(t), loc='right')
    plt.axis('equal')
    plt.draw()
    input("Press Enter to continue...")
    C.remove()

#------------------------------------------------------------------------------

from scipy import *
from matplotlib import pyplot as mpl
from os import sys

##
# Run like
#  $ python viz-ex-04.py serial|braid
#
#  Vizualize the files from ex-04.c or ex-04-serial.c
#       - ex-04.out.state      [ output from ex-04-serial ]
#       - ex-04.out.***.###    [ output from ex-04, *** is the step number, and ### is the proc number ]
#       - ex-04.out.design     [ output by both ex-04 and ex-04-serial ]
#
#  The file ex-04.out.design is a file with Ntime lines, and line k has time
#  step k's design value.
#
#  The file ex-04.out.state is a file with Ntime lines, and line k has time
#  step k's two solution values.  This is only output from ex-04-serial.
#
#  The files ex-04.out.***.###  This file contains one value, the two solution
#  values at the corresponding step number.  This is only output from ex-04.
#
##

# if len(sys.argv) != 2:
#     print " Please run this script with one command line argument, serial|braid, "
#     print " indicating if you want to visualize the output from ex-04-serial or ex-04."
#     sys.exit(1)
#
# ##
# # Load design
# design = loadtxt('ex-04.out.design')
# nsteps = design.shape[0]
# tmesh = linspace(0,1.0,nsteps+1)
#
# ##
# # Load State
# if sys.argv[1] == "serial":
#     state_vec = loadtxt('ex-04.out.state')
# elif sys.argv[1] == "braid":
#     current_rank = 0
#     file_stem = "ex-04.out."
#     state_vec = zeros((nsteps,2))
#     # We don't know the MPI ranks ahead of time that generated the files, so we guess :-)
#     for step in range(1,nsteps+1):
#         try:
#             fname = file_stem + "%04d"%step + '.' + "%03d"%current_rank
#             state_vec[step-1,:] = loadtxt(fname, delimiter=',')
#         except:
#             current_rank = current_rank + 1
#             fname = file_stem + "%04d"%step + '.' + "%03d"%current_rank
#             state_vec[step-1,:] = loadtxt(fname, delimiter=',')
#
# else:
#     print " Please run this script with one command line argument, serial|braid, "
#     print " indicating if you want to visualize the output from ex-04-serial or ex-04."
#     sys.exit(1)
#
# mpl.figure(1)
# mpl.plot(tmesh[1:], state_vec[:,0], '-b')
# mpl.plot(tmesh[1:], state_vec[:,1], '-k')
# mpl.ylabel('u')
# mpl.xlabel('time')
# mpl.title('Solution Values')
# mpl.legend(['Component 1', 'Component 2'])
#
# mpl.figure(2)
# mpl.plot(tmesh[1:], design, '-k')
# mpl.ylabel('design')
# mpl.xlabel('time')
# mpl.title('Design Solution')
# mpl.show()


with open('trischur-adv-diff.out.u.000') as f:
        lines = f.readlines()
for line in lines:
    line = line[6:]
    split = line.split(',')
mspace=len(split)
nsteps=len(lines)
#create the t mesh
tmesh = linspace(0,1.0,nsteps)
xmesh = linspace(0,1.0,mspace)
##
current_rank = 0
state_vec = empty([nsteps, mspace])

with open('trischur-adv-diff.out.u.000') as f:
    lines = f.readlines()
count = 0
for line in lines:
    line = line[6:]
    split = line.split(',')
    count2 = 0
    for thing in split:
        state_vec[count,count2] = float(split[count2])
        count2+=1
    count+=1

# mpl.figure(1)
# for x in range(0,mspace-1):
#     mpl.plot(tmesh[1:], state_vec[:,x], '-b')
# mpl.ylabel('U')
# mpl.xlabel('time')
# mpl.title('Solution Values')

# mpl.figure(2)
# for x in range(0,mspace-1):
#     mpl.plot(tmesh[1:], w_vec[:,x], '-b')
# mpl.ylabel('W')
# mpl.xlabel('time')
# mpl.title('Solution Values')

# mpl.figure(3)
# for x in range(0,mspace-1):
#     mpl.plot(tmesh[1:], v_vec[:,x], '-b')
# mpl.ylabel('V')
# mpl.xlabel('time')
# mpl.title('Solution Values')
# mpl.show()

import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set up grid and test data
nx, ny = mspace, nsteps
x = range(nx)
y = range(ny)


hf = plt.figure(1)
ha = hf.add_subplot(111, projection='3d')

X, Y = numpy.meshgrid(xmesh, tmesh)  # `plot_surface` expects `x` and `y` data to be 2D
ha.plot_surface(X, Y, state_vec)
ha.set_xlabel('X (position)')
ha.set_ylabel('t (time)')
ha.set_zlabel('U value')


plt.show()

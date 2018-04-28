import numpy as np
import h5py
from post_process import PPCLASS
import pdb
import sys
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

"""
Test to see if PPCLASS post_process reads data correctly by recreating mass density 
animation
 
Tests:
1. read and plot mass density for t = 0
2. read and plot mass density for t 
"""

#post process object
pp = PPCLASS()

pp.readMasterFile('/home/samserra/Projects/QPO-hydrodynamic-simulation-analysis/Data/torusrc8p3_0p2eta1_4b3_1')


nStart = 0
nStop = int(pp.mDumpIDArray[-1])
iStepSize = 1

# read mass density for t = 0
time = pp.getCurrentTime(0)
pp.readData(0)

# plot 
fig, ax = plt.subplots()
ax.scatter(pp.mZoneCenter[0::2],pp.mZoneCenter[1::2],c=pp.mScalarFields[3,:],s=1)
ax.set_xlabel('R')
ax.set_ylabel('Z')
ax.set_title('mass density t = {}'.format(pp.mTimeArray[0]))
fig.savefig('Output/massDensityTest.png')

# animate function

def animate(fnum):
    
    time = pp.getCurrentTime(fnum)

    # progress bar
    sys.stdout.write('\rTime {}/{}'.format(pp.mTimeArray[fnum],pp.mTimeArray[-1]))
    sys.stdout.flush()
    
    # read data
    pp.readData(fnum)  #t = 0 is dumpID of 0 (zeroth time step)

    # clear axes
    ax.clear()

    # replot axis title and labels
    ax.set_xlabel('R')
    ax.set_ylabel('Z')
    ax.set_title('mass density t = {}'.format(pp.mTimeArray[fnum]))

    # plot mass density for t = fnum
    ax.scatter(pp.mZoneCenter[0::2],pp.mZoneCenter[1::2],c=pp.mScalarFields[3,:],s=1)

# animate
anim = FuncAnimation(fig, animate, interval = 100, frames = len(pp.mTimeArray))
anim.save('Output/massDensityAnimation.mp4')

print("Got to end")

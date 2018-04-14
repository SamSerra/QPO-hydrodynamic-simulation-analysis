import numpy as np
import h5py
from post_process import PPCLASS
import pdb
import matplotlib.pyplot as plt

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
ax.scatter(pp.mZoneCenter[0::2],pp.mZoneCenter[1::2])
fig.savefig('gridPlot.png')

#ax.contourf(pp.massDensity, cmap="PRGn")
#fig.savefig('massDensityTest.png')




'''
for n in np.arange(nStart,nStop,iStepSize):
    time = pp.getCurrentTime(n)
    print("Time: {} \t Timestep: {}".format(time,n))

    pp.readData(n)  #t = 0 is dumpID of 0 (zeroth time step)
'''

print("Got to end")

import numpy as np
import h5py
from post_process import PPCLASS
import pdb

"""
Test to see if PPCLASS post_process file works 
Tests:
1. readMasterFile works properly ---Done
2. attempt to read velocity data for t = 0 --Done
3. attempt to read velocity data for t --Done

"""

#post process object
pp = PPCLASS()

pp.readMasterFile('/home/samserra/Projects/QPO-hydrodynamic-simulation-analysis/Data/torusrc8p3_0p2eta1_4b3_1')


nStart = 0
nStop = int(pp.mDumpIDArray[-1])
iStepSize = 1

for n in np.arange(nStart,nStop,iStepSize):
    time = pp.getCurrentTime(n)
    print("Time: {} \t Timestep: {}".format(time,n))

    pp.readData(n)  #t = 0 is dumpID of 0 (zeroth time step)

print(pp.velocityField)
# vel data for t=0
pp.readData(0)  #t = 0 is dumpID of 0 (zeroth time step)

print("Got to the end")

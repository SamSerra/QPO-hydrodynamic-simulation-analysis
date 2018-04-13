import numpy as np
import h5py
from post_process import PPCLASS
import pdb

"""
Test to see if PPCLASS post_process file works 
Tests:
1. readMasterFile works properly ---Done
<<<<<<< HEAD
2. attempt to read velocity data for t = 0 --Done
3. attempt to read velocity data for t --Done
=======
2. attempt to read velocity data for t = 0
3. attempt to read velocity data for t 
>>>>>>> cb53c668715043a9b42d341e46eddaaa11e4ce48

"""

#post process object
pp = PPCLASS()

pp.readMasterFile('/home/samserra/Projects/QPO-hydrodynamic-simulation-analysis/Data/torusrc8p3_0p2eta1_4b3_1')

<<<<<<< HEAD

nStart = 0
nStop = int(pp.mDumpIDArray[-1])
iStepSize = 1


for n in np.arange(nStart,nStop,iStepSize):
    time = pp.getCurrentTime(n)
    print("Time: {} \t Timestep: {}".format(time,n))

    pp.readData(n)  #t = 0 is dumpID of 0 (zeroth time step)

print(pp.velocityField)
=======
# read vel data for t=0
pp.readData(0)  #t = 0 is dumpID of 0 (zeroth time step)
print(pp.velocityField)
#print("xdata\n",pp.velocityX)
#print("ydata\n",pp.velocityY)
>>>>>>> cb53c668715043a9b42d341e46eddaaa11e4ce48

print("Got to the end")

import numpy as np
import h5py
from post_process import PPCLASS
import pdb

"""
Test to see if PPCLASS post_process file works 
Tests:
1. readMasterFile works properly ---Done
2. attempt to read velocity data for t = 0
3. attempt to read velocity data for t 

"""

#post process object
pp = PPCLASS()

pp.readMasterFile('/home/samserra/Projects/QPO-hydrodynamic-simulation-analysis/Data/torusrc8p3_0p2eta1_4b3_1')

# read vel data for t=0
pp.readData(0)  #t = 0 is dumpID of 0 (zeroth time step)
print(pp.velocityField)
#print("xdata\n",pp.velocityX)
#print("ydata\n",pp.velocityY)

print("Got to the end")

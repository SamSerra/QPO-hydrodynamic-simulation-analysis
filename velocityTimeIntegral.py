import numpy as np
import h5py
from post_process import PPCLASS
import pdb
import sys 
"""
Integrate x and y velocity data over time at particular frequencies using a trapazoidal scheme.
"""
# Set up simulation parameters
#-----------------------------------------------------------------------------#


# constants (cgs)
solarMass = 1.989e+33 
cee = 2.9979e+10
bigGee = 6.6732e-8

# parameters
bhMass = 10 * solarMass
timeUnit = bigGee * bhMass / cee**3


# Read data
#-----------------------------------------------------------------------------#

#post process object
pp = PPCLASS()

pp.readMasterFile('/home/samserra/Projects/QPO-hydrodynamic-simulation-analysis/Data/torusrc8p3_0p2eta1_4b3_1')


# start and stop times
nStart = 0
nStop = 10#int(pp.mDumpIDArray[-1])
iStepSize = 1

# blank x,y velocity array
xVelArray = np.array([])
yVelArray = np.array([])

# convert time to cgs
timeArray = pp.mTimeArray/timeUnit

# difference array for time integral
dt = np.diff(timeArray)

# read vel data for all times 
for t in np.arange(nStart,nStop,iStepSize):
    time = pp.getCurrentTime(t)
    # progress bar
    sys.stdout.write('\rReading Velocity for Time {}/{}'.format(time, pp.mTimeArray[-1]))
    sys.stdout.flush()

    # read data and get x,y velocities
    pp.readData(t) 

    xVelCurrent, yVelCurrent = pp.velocityField[::2], pp.velocityField[1::2]
    xVelArray, yVelArray = np.append(xVelArray,xVelCurrent, axis=0), np.append(yVelArray,yVelCurrent, axis=0)
    print(xVelArray.shape)
# Integration
#------------------------------------------------------------------------------------#
print('Taking integral...')

# dictionaries for differnt integrals
functionDict = {0:np.cos, 1:np.sin} # cosine and sine
freqDict = {'breathing':203, 'plus':106} # breathing and plus-mode frequencies

# blank array to hold integrals
xVelIntegral = np.array([])
yVelIntegral = np.array([])

# take integral for x and y velocity for 4 different combintations of function and frequency:
# cosine and breathing, cosine and plus, sine and breathing, sine and plus
# total of 8 integrals
for fun in np.arange(len(functionDict)):
    for freq in list(freqDict.keys()):
        
        print('{} and {}'.format(functionDict[fun], freq))
        print(xVelArray.T.shape, timeArray.shape)
        # find heights of trapazoids
        newXVelArray = xVelArray.T * functionDict[fun](2*np.pi*freqDict[freq]*timeArray)
        newYVelArray = yVelArray.T * functionDict[fun](2*np.pi*freqDict[freq]*timeArray)

        # integration scheme
        areaX = np.sum(.5*(newXVelArray[:, ::-1] + newXVelArray[:, 1::]) * dt, axis=1)
        areaY = np.sum(.5*(newYVelArray[::-1,:] + newYVelArray[1::,:]) * dt, axis=1)

        xVelIntegral = np.append(xVelIntegral, areaX)
        yVelIntegral = np.append(yVelIntegral, areaY)

# Make datafile
#-----------------------------------------------------------------------------------#

# fetch grid
R, Z = pp.mZoneCenter[::2], pp.mZoneCenter[1::2]

# write to txt file
with open('Output/velocityTimeIntegral.txt') as f:
    # file header
    f.write('{:>10} | {:>10} | {:>10} | {:>10} | {:>10} | {:>10} | {:>10} | {:>10} | {:>10} | {:>10}'.format('#R', 'Z', 'Vx Cosine B.', 'Vx Cosine +.','Vx Sine B.', 'Vx Sine +.','Vy Cosine B.', 'Vy Cosine +.','Vy Sine B.', 'Vy Sine +.'))

    # write data
    for i in np.arange(len(R)):
        f.write('{:>10} | {:>10} | {:>10} | {:>10} | {:>10} | {:>10} | {:>10} | {:>10} | {:>10} | {:>10}'.format(R[i], Z[i], xVelIntegral[0,i], xVelIntegral[1,i], xVelIntegral[2,i], xVelIntegral[3,i], yVelIntegral[0,i], yVelIntegral[1,i], yVelIntegral[2,i], yVelIntegral[3,i]))


print("Got to the end")

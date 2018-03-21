import numpy as np
import matplotlib.pyplot as plt
"""
Module to Supplement main_pp.py
"""

def ConvertToPolar(nDims, xArray, yArray, zArray = np.zeros(1,)):
        """
        Takes in arrays of cartesian coordinates and transforms
        them to polar or spherical coordinates.
        Spherical coordinates set up as rho, phi, and theta (azimuth)

        Returns radius array, phi, theta array
        """
        assert xArray.shape == yArray.shape, "Array's must have same shape!"

        if nDims == 3:
                pass
        elif nDims == 2:
                zArray=np.zeros(xArray.shape)
                
        #create arrays for coordinates
        rho, phi, theta = np.ones(xArray.shape), np.ones(xArray.shape), np.ones(xArray.shape)
        
        rho=np.sqrt(xArray**2+yArray**2+zArray**2) #if 2-dimensional, rho == rad
        rad=np.sqrt(xArray**2+yArray**2)
        
        phi=np.arctan2(rad,zArray)
        
        theta=np.arctan2(yArray,xArray)

        return rho, phi, theta

def plot_polar_contour(values, theta, radius):
        """
        Takes in z,theta,r arrays
        Generates a polar contour plot symetric about theta = 0
        """
        theta = np.array(theta)
        radius = np.array(radius)

        values = np.array(values)
        values = values.reshape(len(theta), len(radius))

        r, tet = np.meshgrid(radius, theta)	#meshgrid set up opposite of values b/c pyplot uses transpose 

        print("Generating Plot...")	
        maxbound = np.max(theta)	#bounds for plot

        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        cax = ax.contourf(tet, r, values, 30, cmap = "hot", vmin = 1.0e-29, vmax = 1.0e-18)

        cb = fig.colorbar(cax)
        cb.set_label("Mass Density")

        x=ax.axes.get_xaxis()
        x.set_visible(False)

        plt.savefig('MassDensityContour.png', dpi = 600)
        print("Plot Saved")
        return fig, ax, cax

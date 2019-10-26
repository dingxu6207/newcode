# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 02:32:53 2019

@author: dingxu
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture


hdu = fits.open('YFCf050164.fits')
imgdata = hdu[1].data  #hdu[0].header

hdu1 = fits.open('YFCf050165.fits')
imgdata1 = hdu1[1].data  #hdu[0].header

def adjustimage(imagedata, coffe=3):
    mean = np.mean(imagedata)
    sigma = np.std(imagedata)
    mindata = np.min(imagedata)
    maxdata = np.max(imagedata)
    Imin = mean - coffe*sigma
    Imax = mean + coffe*sigma
        
    mindata = max(Imin,mindata)
    maxdata = min(Imax,maxdata)
    return mindata,maxdata

mindata,maxdata = adjustimage(imgdata)

def findsource(imgdata):
    mean, median, std = sigma_clipped_stats(imgdata, sigma=3.0) 
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
    sources = daofind(imgdata - median)
    return sources

sources = findsource(imgdata)
sources1 = findsource(imgdata1)

def printsources(sources):
    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output
    print(sources)




positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
positions1 = np.transpose((sources1['xcentroid'], sources1['ycentroid']))
apertures = CircularAperture(positions, r=10.)

plt.imshow(imgdata1, cmap='gray', vmin = mindata, vmax = maxdata, origin='lower')
plt.plot(647.48, 633.23,'*' )   #647.765	632.988

#apertures.plot(color='blue', lw=1.5, alpha=0.5)
#plt.imshow(imgdata, vmin = mindata, vmax = maxdata, cmap = 'gray')

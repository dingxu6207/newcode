# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 11:02:43 2019

@author: dingxu
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture


fitsname1 = 'E:\\AST3\\RA0252DEC3212\\'+'L20161102_04746_025202+3212_60S_SI_182211.FITS'
fitsname2 = 'E:\\AST3\\RA0252DEC3212\\'+'L20180310_08663_025202+3212_60S_SI_278406.FITS'

onehdu = fits.open(fitsname1)
imgdata1 = onehdu[0].data  #hdu[0].header

twohdu = fits.open(fitsname2)
imgdata2 = twohdu[0].data  #hdu[0].header

def adjustimage(imagedata, coffe):
    mean = np.mean(imagedata)
    sigma = np.std(imagedata)
    mindata = np.min(imagedata)
    maxdata = np.max(imagedata)
    Imin = mean - coffe*sigma
    Imax = mean + coffe*sigma
        
    mindata = max(Imin,mindata)
    maxdata = min(Imax,maxdata)
    return mindata,maxdata


def displayimage(img, coff, i):
    minimg,maximg = adjustimage(img, coff)
    plt.figure(i)
    plt.imshow(img, cmap='gray', vmin = minimg, vmax = maximg)
    plt.plot(705.918, 510.251,'*')

def findsource(img):    
    mean, median, std = sigma_clipped_stats(img, sigma=3.0) 
    daofind = DAOStarFinder(fwhm=5.0, threshold=4.*std)
    sources = daofind(img - median)
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    pmag = np.transpose((sources['xcentroid'], sources['ycentroid'],sources['peak'], sources['mag'],sources['flux']))
    return pmag,positions

pmag1, positions1 = findsource(imgdata1)
displayimage(imgdata1,3,0)
apertures1 = CircularAperture(positions1, r=10.)
apertures1.plot(color='blue', lw=1.5, alpha=0.5)
#print(pmag1)

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 11:02:43 2019

@author: dingxu
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

fitsname1 = 'E:\\BOOTES4\\20181118\\03095\\'+'20181118130453-727-RA.fits'
fitsname2 = 'E:\\BOOTES4\\20181118\\03095\\'+'20181118125001-285-RA.fits'


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


displayimage(imgdata1,3,0)
displayimage(imgdata2,3,1)

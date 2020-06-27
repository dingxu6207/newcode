# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 16:59:45 2020

@author: dingxu
"""

import astroalign as aa
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture
import cv2
import os


fitsname1 = 'E:/AST3/xuyi/'+'L20190406_10053_202404 0412_60S_VR_335648.FITS'
fitsname2 = 'E:/AST3/xuyi/'+'L20190406_10053_202408 0411_60S_VR_335680.FITS'


onehdu = fits.open(fitsname1)
imgdata1 = onehdu[0].data  #hdu[0].header
oneimgdata = np.copy(imgdata1)

twohdu = fits.open(fitsname2)
imgdata2 = twohdu[0].data  #hdu[0].header
twoimgdata = np.copy(imgdata2)

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
    plt.savefig(str(i)+'.jpg')
    

img_aligned, footprint = aa.register(oneimgdata, twoimgdata)
displayimage(oneimgdata,3,1)
displayimage(twoimgdata,3,2)
displayimage(img_aligned,3,3)
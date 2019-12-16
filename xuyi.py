# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 11:27:46 2019

@author: dingxu
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture
import scipy.signal as signal

fitshdu = fits.open('E:\\BOOTES4\\20190606\\'+'201906061524450716.fit')
fitsdata = fitshdu[0].data
fitsdata = np.float32(fitsdata)
fitsdata = signal.medfilt2d(fitsdata, kernel_size=9)  # 二维中值滤波

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



def findsource(img):    
    mean, median, std = sigma_clipped_stats(img, sigma=3.0)
    daofind = DAOStarFinder(fwhm=5, threshold=5.*std)
    sources = daofind(img - median)

    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output
        #print(sources)

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    positionflux = np.transpose((sources['xcentroid'], sources['ycentroid'],  sources['flux']))
    mylist = positionflux.tolist()
    
    return sources,positions,mylist

sources1,positions1,mylist1 =  findsource(fitsdata)
apertures1 = CircularAperture(positions1, r=10.)

mindata,maxdata = adjustimage(fitsdata, 3)
plt.imshow(fitsdata, cmap='gray', vmin = mindata, vmax = maxdata, origin='lower')
apertures1.plot(color='blue', lw=1.5, alpha=0.5)
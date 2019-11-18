# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 14:10:43 2019

@author: dingxu
"""

from astropy.io import fits
import sim
import numpy as np
import matplotlib.pyplot as plt

fitsname1 = 'd8549.0091.fits'
fitsname2 = 'd8549.0088.fits'
routename1 = 'E:\\shunbianyuan\\code\\phot\\'
routename2 = 'E:\\shunbianyuan\\code\\phot\\'

fitsname1 = routename1+fitsname1
fitsname2 = routename2+fitsname2


onehdu = fits.open(fitsname1)
oneimgdata = onehdu[1].data  #hdu[0].header

twohdu = fits.open(fitsname2)
twoimgdata = twohdu[1].data  #hdu[0].header


dx0,dy0=0,0
for i in range(1):
    dx,dy,cor=sim.xcorrcenter(oneimgdata,twoimgdata)
    print(i,dx,dy)
    twoimgdata=sim.immove(twoimgdata,dx,dy)
    dx0=np.max((dx,dx0))
    dy0=np.max((dy,dy0))
    
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


mindata1,maxdata1 = adjustimage(oneimgdata,1)
mindata2,maxdata2 = adjustimage(twoimgdata,1)
    
plt.figure(0)
plt.imshow(oneimgdata, cmap='gray', vmin = mindata1, vmax = maxdata1, origin='lower')

plt.figure(1)
plt.imshow(twoimgdata, cmap='gray', vmin = mindata2, vmax = maxdata2, origin='lower')

plt.figure(3)  
subimg = np.float64(oneimgdata) - np.float64(twoimgdata)
minsub,maxsub = adjustimage(subimg,1)
plt.imshow(subimg, cmap='gray', vmin = minsub, vmax = maxsub, origin='lower')

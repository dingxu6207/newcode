# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 02:58:54 2019

@author: dingxu
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import create_matching_kernel
from photutils import TopHatWindow


hdu = fits.open('YFCf050162.fits')
imgdata = hdu[1].data  #hdu[0].header

hdu1 = fits.open('YFCf050169.fits')
imgdata1 = hdu1[1].data  #hdu[0].header

Cutoutimg = imgdata
Cutoutimg1 = imgdata1

def make_blurred(gray, PSF, eps=0.000):
    input_fft = np.fft.fft2(gray)# 进行二维数组的傅里叶变换
    PSF_fft = np.fft.fft2(PSF)+ eps
    blurred = np.fft.ifft2(input_fft * PSF_fft)
    blurredimage = np.abs(np.fft.fftshift(blurred))
    
    return blurredimage


def adjustimage(imagedata, coffe=4):
    mean = np.mean(imagedata)
    sigma = np.std(imagedata)
    mindata = np.min(imagedata)
    maxdata = np.max(imagedata)
    Imin = mean - coffe*sigma
    Imax = mean + coffe*sigma
        
    mindata = max(Imin,mindata)
    maxdata = min(Imax,maxdata)
    return mindata,maxdata

mindata,maxdata = adjustimage(Cutoutimg)
mindata1,maxdata1 = adjustimage(Cutoutimg1)

plt.figure(0)
plt.imshow(Cutoutimg,vmin = mindata, vmax = maxdata,cmap = 'gray' )

plt.figure(1)
plt.imshow(Cutoutimg1, vmin = mindata1, vmax = maxdata1,cmap = 'gray' )

plt.figure(2)
doubleimg = np.float64(Cutoutimg)-np.float64(Cutoutimg1)
doubleimg = np.abs(doubleimg)
mindata2,maxdata2 = adjustimage(doubleimg)
plt.imshow(doubleimg,vmin=mindata2,vmax=maxdata2, cmap = 'gray' )

window = TopHatWindow(0.8)
fimg = np.float64(Cutoutimg)
fimg1 = np.float64(Cutoutimg1)
kernel = create_matching_kernel(fimg, fimg1, window)

bluimage = make_blurred(Cutoutimg,kernel)
plt.figure(3)
mindata3,maxdata3 = adjustimage(bluimage)
plt.imshow(bluimage,vmin=mindata3,vmax=maxdata3,cmap = 'gray' )

plt.figure(4)
doubleimg1 = np.float64(Cutoutimg1)-np.float64(bluimage)
doubleimg1 = np.abs(doubleimg1)
mindata4,maxdata4 = adjustimage(doubleimg1)
plt.imshow(doubleimg1,vmin=mindata4,vmax=maxdata4,cmap = 'gray' )

# -*- coding: utf-8 -*-
"""
Created on Sun May 10 14:59:54 2020

@author: dingxu
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture
from sklearn.decomposition import PCA
import cv2
import os
import math

fitsname1 = 'E:\\shunbianyuan\\newdata\\'+'201906061614530716.fit'
fitsname2 = 'E:\\shunbianyuan\\newdata\\'+'201906061616150716.fit'
onehdu = fits.open(fitsname1)
imgdata1 = onehdu[0].data  #hdu[0].header
copydata1 = np.copy(imgdata1)
imgdata1 = np.float32(copydata1)
oneimgdata = imgdata1
timedate = onehdu[0].header['DATE']
#oneimgdata = signal.medfilt2d(imgdata1, kernel_size=5)  # 二维中值滤波
hang1,lie1 = oneimgdata.shape

twohdu = fits.open(fitsname2)
imgdata2 = twohdu[0].data  #hdu[0].header
copydata2 = np.copy(imgdata2)
imgdata2 = np.float32(copydata2)
twoimgdata = imgdata2
#twoimgdata = signal.medfilt2d(imgdata2, kernel_size=5)  # 二维中值滤波
hang2,lie2 = twoimgdata.shape


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


def findsource(img):    
    mean, median, std = sigma_clipped_stats(img, sigma=3.0) 
    daofind = DAOStarFinder(fwhm=14, threshold=6.*std)
    sources = daofind(img - median)

    #tezhen = np.transpose((sources['xcentroid'], sources['ycentroid']))
    tezhen = np.transpose((sources['xcentroid'], sources['ycentroid']))
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

    return tezhen,positions


###实现找星###
tezhen1,positions1 =  findsource(oneimgdata)
tezhen2,positions2 =  findsource(twoimgdata)
   
apertures1 = CircularAperture(positions1, r=7.)
apertures2 = CircularAperture(positions2, r=7.)

displayimage(oneimgdata,3,0)
apertures1.plot(color='blue', lw=1.5, alpha=0.5)

displayimage(twoimgdata,3,1)
apertures2.plot(color='blue', lw=1.5, alpha=0.5)

x1 = int(positions2[3][0]-8)
x2 = int(positions2[3][0]+8) 

y1 = int(positions2[3][1]-8)
y2 = int(positions2[3][1]+8) 
displayimage(twoimgdata[y1:y2,x1:x2],3,2)
flatdata = twoimgdata[y1:y2,x1:x2]

pca = PCA(n_components=5)
pca.fit(flatdata)
pcadata = pca.transform(flatdata)
displayimage(pcadata,3,3)
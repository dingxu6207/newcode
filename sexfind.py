# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 15:51:49 2020

@author: dingxu
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture

filename = 'E:\\shunbianyuan\\newdata\\'+'M31.fts'
filecat = 'E:\\shunbianyuan\\newdata\\test.cat'
fitshdu = fits.open(filename)
fitsdata = fitshdu[0].data

fdata = np.loadtxt(filecat)
hang,lie = fdata.shape
matrixdata = fdata[0:hang, 1:lie]

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
    
displayimage(fitsdata, 3, 0)

for i in range(hang):
    xmin = matrixdata[i,0]-5
    ymin = matrixdata[i,1]-5
    xmax = matrixdata[i,2]+5
    ymax = matrixdata[i,3]+5                   
    plt.plot([xmin,xmax],[ymin,ymin],linewidth = 1.0,color = 'b')
    plt.plot([xmin,xmin],[ymin,ymax],linewidth = 1.0,color = 'b')
    plt.plot([xmin,xmax],[ymax,ymax],linewidth = 1.0,color = 'b')
    plt.plot([xmax,xmax],[ymin,ymax],linewidth = 1.0,color = 'b')

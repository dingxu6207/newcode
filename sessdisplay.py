# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 22:58:01 2020

@author: dingxu
"""

import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry

file  = '0.fits'
path = 'E:\\shunbianyuan\\dataxingtuan\\alngc7142\\'
filename = path+file
fitshdu = fits.open(filename)
data = fitshdu[0].data
i = 2 #行扫描 i = 21
j = 3 #列扫描 j=20
fitsdata = data[398*i:398+398*i,389*j:389+389*j]

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
    
    

starlight = np.loadtxt('starlight.txt')
hang,lie = starlight.shape

targetx = 326
targety = 134
threshold = 10
for i in range(hang):
    delt = np.sqrt((targetx - starlight[i][0])**2+(targety - starlight[i][1])**2)
    if delt < threshold:
        print(i)
        m = i
          
x = np.arange(0,lie-2,1,dtype=np.int16)

for i in range(hang):
    plt.clf()
    displayimage(fitsdata, 1, 0)
    plt.plot(starlight[i,0], starlight[i,1],'*')
    plt.pause(3)
    
    plt.clf()
    plt.figure(1)        
    fluxy = starlight[i,2:]-starlight[m,2:]
    plt.plot(x,fluxy,'.')
    
    nihe = np.polyfit(x, fluxy, 0)
    nihey = np.poly1d(nihe)
    y = nihey(x)
    plt.plot(x, y, label="xx")
    
    plt.hlines(np.mean(fluxy)+2*np.std(fluxy), 0, lie,color="red")#横线
    plt.hlines(np.mean(fluxy)-2*np.std(fluxy), 0, lie,color="red")#横线
   
    plt.legend()
    plt.title(str(i))
    plt.pause(3)
    


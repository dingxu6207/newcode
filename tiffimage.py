# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 20:54:28 2020

@author: dingxu
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import imageio

filetemp = []
count = 0
oripath = 'E:\\shunbianyuan\\todingx\\aligendata\\'
for root, dirs, files in os.walk(oripath):
   for file in files:
       #print(file)
       if (file[-1] == 's'):
           count = count+1
           filetemp.append(file)
       
       
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
    #plt.plot(376.82, 137.232, '*')
    #plt.savefig(oripath+str(i)+'.jpg')

i = 0 
for i in range(0, count):
    try:
        fitshdu = fits.open(oripath+filetemp[i])
        fitsdata = fitshdu[0].data 
        displayimage(fitsdata, 1,1)
        plt.pause(0.001)
    except:
        print('error')
  
'''
gif_images = []   
i = 0   
for i in range(0, count):
    try:
        print(oripath+str(i)+'.jpg')
        gif_images.append(imageio.imread(oripath+str(i)+'.jpg'))
    except:
        print('eroor!!!')
        
imageio.mimsave(oripath+"test.gif",gif_images,fps=1)        
'''
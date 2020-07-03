# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 15:37:54 2020

@author: dingxu
"""

import os
import astroalign as aa
import numpy as np
from astropy.io import fits

filetemp = []
count = 0
oripath = 'E:\\shunbianyuan\\todingx\\origindata\\'
for root, dirs, files in os.walk(oripath):
   for file in files:
       #print(file)
       count = count+1
       filetemp.append(file)
       
def witefits(data,name):
    writepath = 'E:\\shunbianyuan\\todingx\\aligendata\\'
    name = str(name)
    os.chdir(writepath)
    if(os.path.exists(writepath+name + '.fits')):
       os.remove(name + '.fits')
    grey=fits.PrimaryHDU(data)
    greyHDU=fits.HDUList([grey])
    greyHDU.writeto(name + '.fits')
 

fitsname1 = oripath+filetemp[0]
onehdu = fits.open(fitsname1)
imgdata1 = onehdu[0].data  #hdu[0].header
copydata1 = np.copy(imgdata1)   
witefits(copydata1,0)  
 
for i in range(1, count):
    fitsname2 = oripath+filetemp[i]
    twohdu = fits.open(fitsname2)
    imgdata2 = twohdu[0].data  #hdu[0].header
    copydata2 = np.copy(imgdata2)  
    try:
        aligned_image, footprint = aa.register(copydata2, copydata1)
        witefits(aligned_image,i)
        print('ok!')
    except:
        print('error!!!')
    
    
       
       
       
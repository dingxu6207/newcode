# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 21:07:44 2020

@author: dingxu
"""

import os
import numpy as np
from astropy.io import fits
import sajiaofunc
#import imageio

filetemp = []
count = 0
oripath = 'E:\\shunbianyuan\\dataxingtuan\\ngc7423\\'
for root, dirs, files in os.walk(oripath):
   for file in files:
       if (file[-4:] == '.fit'):
       #print(file)
           count = count+1
           filetemp.append(file)
   
   
def witefits(data,name): 
    writepath = 'E:\\shunbianyuan\\dataxingtuan\\alngc7142\\'
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
onedata = np.float32(copydata1)  
witefits(copydata1,0)  
 

#temptiff = []
for i in range(1, count):
    fitsname2 = oripath+filetemp[i]
    twohdu = fits.open(fitsname2)
    imgdata2 = twohdu[0].data  #hdu[0].header
    copydata2 = np.copy(imgdata2)  
    twodata = np.float32(copydata2)
    try:
        aligin = sajiaofunc.newdata(twodata, onedata)
        witefits(aligin,i)
        #temptiff.append(aligned_image)
        print('ok!')
    except:
        print('error!!!')
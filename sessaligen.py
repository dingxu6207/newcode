# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 15:37:54 2020

@author: dingxu
"""

import os
import astroalign as aa
import numpy as np
from astropy.io import fits
#import imageio

filetemp = []
count = 0
oripath = 'E:\\shunbianyuan\\dataxingtuan\\ngc7142\\'
for root, dirs, files in os.walk(oripath):
   for file in files:
       if (file[-4:] == '.fit'):
       #print(file)
           count = count+1
           filetemp.append(file)
      
def witefits(data,name,timedate):
    writepath = 'E:\\shunbianyuan\\dataxingtuan\\alngc7142\\'
    name = str(name)
    os.chdir(writepath)
    if(os.path.exists(writepath+name + '.fits')):
        os.remove(name + '.fits')
    hdr = fits.Header()
    hdr['DATE'] = timedate
    fitsdata = np.float32(data)
    fits.writeto(name + '.fits', fitsdata, hdr)

fitsname1 = oripath+filetemp[1]
onehdu = fits.open(fitsname1)
datatime = onehdu[0].header['DATE']
imgdata1 = onehdu[0].data  #hdu[0].header
copydata1 = np.copy(imgdata1)   
witefits(copydata1,0,datatime)  
 

for i in range(1, count):
    fitsname2 = oripath+filetemp[i]
    twohdu = fits.open(fitsname2)
    imgdata2 = twohdu[0].data  #hdu[0].header
    datatime2 = twohdu[0].header['DATE']
    #print(datatime2,fitsname2)
    copydata2 = np.copy(imgdata2)  
    try:
        aligned_image, footprint = aa.register(copydata2, copydata1)
        witefits(aligned_image, i, datatime2)
        print('ok!')
    except:
        print('error!!!')
    


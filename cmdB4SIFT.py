# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 15:40:52 2019

@author: dingxu
"""

import os
from astropy.io import fits
import numpy as np

curentpath = 'E:\\BOOTES4\\20190606\\'
filelist = []
for root, dirs, files in os.walk(curentpath):
   for file in files:
       strfile = os.path.join(root, file)
       if (strfile[-5:] == '6.fit'):
           hdu = fits.open(strfile)
           data = hdu[0].data
           meandata = np.mean(data)
           #print(strfile,meandata)
           filelist.append(file)



length = len(filelist)
name2 = filelist[0]

for i in range(1,length):
    name1 = filelist[i]
    cmd = 'python'+' '+'B4sift.py'+' '+name1+' '+name2+' '+ str(i)
    os.system(cmd)

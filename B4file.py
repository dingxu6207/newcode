a# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 09:34:59 2019

@author: dingxu
"""

import os
from astropy.io import fits
import numpy as np
import shutil

os.chdir('J:\\bootes4\\xdr\\BOOTES-4')
curentpath = os.getcwd()
print(curentpath)

dirpath = 'J:\\0716'

i = 0
for root, dirs, files in os.walk(curentpath):
   for file in files:
       strfile = os.path.join(root, file)
       if (strfile[-5:] == '6.fit'):
           hdu = fits.open(strfile)
           data = hdu[0].data
           meandata = np.mean(data)
           print(strfile,meandata)
           shutil.copy(strfile,dirpath)
           i = i+1
print(i)
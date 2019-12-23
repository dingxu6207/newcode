# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 09:34:59 2019

@author: dingxu
"""

import os
from astropy.io import fits
import numpy as np

os.chdir('E:\\BOOTES4\\20181118\\03095')
curentpath = os.getcwd()
print(curentpath)



for root, dirs, files in os.walk(curentpath):
   for file in files:
       strfile = os.path.join(root, file)
       if (strfile[-5:] == '.fits'):
           hdu = fits.open(strfile)
           data = hdu[0].data
           meandata = np.mean(data)
           print(strfile,meandata)
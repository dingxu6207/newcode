# -*- coding: utf-8 -*-
"""
Created on Sat May  9 17:08:46 2020

@author: dingxu
"""

import matplotlib.pyplot as plt
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

fitsname1 = 'E:\\shunbianyuan\\newdata\\'+'3.fits'

w = WCS('fitsname1')  

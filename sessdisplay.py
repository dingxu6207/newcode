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


starlight = np.loadtxt('starlight.txt')
hang,lie = starlight.shape

targetx = 250
targety = 384
threshold = 10
for i in range(hang):
    delt = np.sqrt((targetx - starlight[i][0])**2+(targety - starlight[i][1])**2)
    if delt < threshold:
        print(i)
            

for i in range(hang):
    plt.figure(1)
    plt.plot(starlight[i,2:]-starlight[55,2:])
    plt.pause(0.01)


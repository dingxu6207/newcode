# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 00:33:52 2019

@author: dingxu
"""
import numpy as np
from astropy.modeling.models import Gaussian2D
import matplotlib.pyplot as plt
from photutils import create_matching_kernel
from photutils import TopHatWindow

y, x = np.mgrid[0:51, 0:51]
gm1 = Gaussian2D(100, 25, 25, 3, 3)
gm2 = Gaussian2D(100, 25, 25, 6, 6)
gm3 = Gaussian2D(120, 25, 25, 5, 5)


def make_blurred(gray, PSF, eps=0):
    input_fft = np.fft.fft2(gray)# 进行二维数组的傅里叶变换
    PSF_fft = np.fft.fft2(PSF)+ eps
    blurred = np.fft.ifft2(input_fft * PSF_fft)
    blurredimage = np.abs(np.fft.fftshift(blurred))
    
    return blurredimage

g1 = gm1(x, y)
g2 = gm2(x, y)
g3 = gm3(x, y)

g1 /= g1.sum()
g2 /= g2.sum()
g3 /= g3.sum()


window = TopHatWindow(0.3)
kernel = create_matching_kernel(g1, g2, window)

bluimage = make_blurred(g1,kernel)

plt.figure(0)
plt.imshow(g1, cmap='gray', origin='lower')

plt.figure(1)
plt.imshow(g2, cmap='gray', origin='lower')

gisubg2 = np.float64(g1) - np.float64(g2)
plt.figure(2)
plt.imshow(gisubg2, cmap='gray', origin='lower')

plt.figure(3)
plt.imshow(bluimage, cmap='gray', origin='lower')

subimage = np.float64(g3) - np.float64(bluimage)
plt.figure(4)
plt.imshow(subimage, cmap='gray', origin='lower')

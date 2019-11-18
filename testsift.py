# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 20:00:23 2019

@author: dingxu
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 15:20:11 2019

@author: dingxu
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture
import cv2



fitsname1 = 'd8549.0088.fits'
fitsname2 = 'd8549.0091.fits'

onehdu = fits.open(fitsname1)
oneimgdata = onehdu[1].data  #hdu[0].header
hang1,lie1 = oneimgdata.shape

twohdu = fits.open(fitsname2)
twoimgdata = twohdu[1].data  #hdu[0].header
hang2,lie2 = twoimgdata.shape



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



def findsource(img):    
    mean, median, std = sigma_clipped_stats(img, sigma=3.0) 
    daofind = DAOStarFinder(fwhm=2.5, threshold=5.*std)
    sources = daofind(img - median)

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

    return positions


###实现找星###
positions1 =  findsource(oneimgdata)
mindata1,maxdata1 = adjustimage(oneimgdata,3)

positions2 =  findsource(twoimgdata)
mindata2,maxdata2 = adjustimage(twoimgdata,3)
    
apertures1 = CircularAperture(positions1, r=10.)
apertures2 = CircularAperture(positions2, r=10.)

plt.figure(0)
plt.imshow(oneimgdata, cmap='gray', vmin = mindata1, vmax = maxdata1)
apertures1.plot(color='blue', lw=1.5, alpha=0.5)

plt.figure(1)
plt.imshow(twoimgdata, cmap='gray', vmin = mindata2, vmax = maxdata2)
plt.savefig('two.jpg')
apertures2.plot(color='blue', lw=1.5, alpha=0.5)


def u8img(img):
    maximg = np.max(img)
    minimg = np.min(img)
    guiyiimg = 255*(img-minimg)/(maximg-minimg)
    return np.uint8(guiyiimg)



uimg1 = u8img(oneimgdata) 
uimg2 = u8img(twoimgdata)



lenposition1 = len(positions1)
lenposition2 = len(positions2)
keyimg1 = np.zeros((lenposition1,128),dtype = np.float32)
keyimg2 = np.zeros((lenposition2,128),dtype = np.float32)
i = 0
j = 0
for i in range(lenposition1):
    keyimg1[i,0:2] = positions1[i,:]
    
for j in range(lenposition2):
    keyimg2[j,0:2] = positions2[j,:]   


# FLANN 参数设计
FLANN_INDEX_KDTREE = 0
index_params = dict(algorithm = FLANN_INDEX_KDTREE, trees = 5)
search_params = dict(checks=50)
flann = cv2.FlannBasedMatcher(index_params,search_params)


matches = flann.knnMatch(keyimg1,keyimg2,k=2)
lenmatches = len(matches)
temp = np.zeros((int(lenmatches/4),2),dtype = np.uint)
j = 0
for i, (m1, m2) in enumerate(matches):
    if m1.distance < 0.4 * m2.distance:# 两个特征向量之间的欧氏距离，越小表明匹配度越高。
        print(i,m1.queryIdx,m1.trainIdx)
        temp[j][0] = m1.queryIdx
        temp[j][1] = m1.trainIdx
        j = j+1
'''
plt.figure(2)
plt.imshow(oneimgdata, cmap='gray', vmin = mindata1, vmax = maxdata1)
plt.plot(positions1[15][0],positions1[15][1],'*')

plt.figure(3)
plt.imshow(twoimgdata, cmap='gray', vmin = mindata2, vmax = maxdata2)
plt.plot(positions2[10][0],positions2[10][1],'*')
'''

hmerge = np.hstack((oneimgdata, twoimgdata)) #水平拼接
minhmerge,maxhmerge = adjustimage(hmerge,3)
plt.figure(2)
plt.imshow(hmerge, cmap='gray', vmin = minhmerge, vmax = maxhmerge)

for i in range(int(lenmatches/4)):
    x = temp[i][0]
    y = temp[i][1]
    x10 = positions1[x][0]
    y10 = positions1[x][1]
    
    x11 = positions2[y][0]
    y11 = positions2[y][1]
    
    plt.plot(x10,y10,'*')
    plt.plot(x11+lie1,y11,'*')
    
    plt.plot([x10,x11+lie1],[y10,y11])


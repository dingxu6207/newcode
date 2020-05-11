# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 21:54:40 2019

@author: dingxu
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture
import cv2
#import scipy.signal as signal
import os
import math

fitsname1 = 'E:\\BOOTES4\\20181118\\03088\\'+'20181118130518-952-RA.fits'
fitsname2 = 'E:\\BOOTES4\\20181118\\03088\\'+'20181118130621-081-RA.fits'

onehdu = fits.open(fitsname1)
imgdata1 = onehdu[0].data  #hdu[0].header
copydata1 = np.copy(imgdata1)
imgdata1 = np.float32(copydata1)
oneimgdata = imgdata1
#oneimgdata = signal.medfilt2d(imgdata1, kernel_size=5)  # 二维中值滤波
hang1,lie1 = oneimgdata.shape

twohdu = fits.open(fitsname2)
imgdata2 = twohdu[0].data  #hdu[0].header
copydata2 = np.copy(imgdata2)
imgdata2 = np.float32(copydata2)
twoimgdata = imgdata2
#twoimgdata = signal.medfilt2d(imgdata2, kernel_size=5)  # 二维中值滤波
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

def displayimage(img, coff, i):
    minimg,maximg = adjustimage(img, coff)
    plt.figure(i)
    plt.imshow(img, cmap='gray', vmin = minimg, vmax = maximg)
    plt.savefig(str(i)+'.jpg')


def findsource(img):    
    mean, median, std = sigma_clipped_stats(img, sigma=3.0) 
    daofind = DAOStarFinder(fwhm=8.5, threshold=5.*std)
    sources = daofind(img - median)

    tezhen = np.transpose((sources['xcentroid'], sources['ycentroid']))
    #tezhen = np.transpose((sources['xcentroid'], sources['ycentroid'],sources['sharpness']))
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

    return tezhen,positions


###实现找星###
tezhen1,positions1 =  findsource(oneimgdata)
tezhen2,positions2 =  findsource(twoimgdata)
   
apertures1 = CircularAperture(positions1, r=5.)
apertures2 = CircularAperture(positions2, r=5.)

displayimage(oneimgdata,3,0)
apertures1.plot(color='blue', lw=1.5, alpha=0.5)

displayimage(twoimgdata,3,1)
apertures2.plot(color='blue', lw=1.5, alpha=0.5)

lenposition1 = len(positions1)
lenposition2 = len(positions2)
keyimg1 = np.zeros((lenposition1,128),dtype = np.float32)
keyimg2 = np.zeros((lenposition2,128),dtype = np.float32)
i = 0
j = 0
for i in range(lenposition1):
    keyimg1[i,0:2] = tezhen1[i,:]
    
for j in range(lenposition2):
    keyimg2[j,0:2] = tezhen2[j,:]   


# FLANN 参数设计
FLANN_INDEX_KDTREE = 0
index_params = dict(algorithm = FLANN_INDEX_KDTREE, trees = 5)
search_params = dict(checks=50)
flann = cv2.FlannBasedMatcher(index_params,search_params)

matches = flann.knnMatch(keyimg1,keyimg2,k=2)

lenpipei = 0
temp1 = []
temp2 = []
for i, (m1, m2) in enumerate(matches):
    if m1.distance < 0.75 * m2.distance:# 两个特征向量之间的欧氏距离，越小表明匹配度越高。
        lenpipei = lenpipei+1
        temp1.append(m1.queryIdx)
        temp2.append(m1.trainIdx)

hmerge = np.hstack((oneimgdata, twoimgdata)) #水平拼接
displayimage(hmerge, 3, 2)

srckp1 = []
srckp2 = []
for i in range(lenpipei):
    x = temp1[i]
    y = temp2[i]
    x10 = positions1[x][0]
    y10 = positions1[x][1]
    srckp1.append(x10)
    srckp1.append(y10)
    src_pts = np.float32(srckp1).reshape(-1,2)
    
    x11 = positions2[y][0]
    y11 = positions2[y][1]
    srckp2.append(x11)
    srckp2.append(y11)
    dst_pts = np.float32(srckp2).reshape(-1,2)
    
    #plt.plot(x10,y10,'*')
    #plt.plot(x11+lie1,y11,'*')
    
    plt.plot([x10,x11+lie1],[y10,y11],linewidth = 0.8)




H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC,5.0)
newimg1 = cv2.warpPerspective(imgdata1, H, (lie1,hang1))

addimg = np.float32(newimg1) + np.float32(imgdata2)
minusimg = np.float32(newimg1) - np.float32(imgdata2)
displayimage(addimg, 3, 3)
displayimage(minusimg, 3, 4)


def witefits(data,name):
    os.remove(name + '.fits')
    grey=fits.PrimaryHDU(data)
    greyHDU=fits.HDUList([grey])
    greyHDU.writeto(name + '.fits')
    
witefits(newimg1,'one')   
witefits(imgdata2,'two') 
witefits(minusimg,'minusimg') 

tempmatrix = np.zeros((3,1),dtype = np.float64)
tempmatrix[2] = 1
deltemp = []

for j in range(lenpipei):
    tempmatrix[0] = src_pts[j][0]
    tempmatrix[1] = src_pts[j][1]
    
    result = np.dot(H,tempmatrix)
    
    rx11 = result[0]/result[2]
    ry11 = result[1]/result[2]
    
    delcha = math.sqrt((rx11-dst_pts[j][0])**2 + (ry11-dst_pts[j][1])**2)
    deltemp.append(delcha)
    
plt.figure(5)
plt.plot(deltemp)
print(np.mean(deltemp[15:40]))    

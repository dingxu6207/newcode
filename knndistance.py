# -*- coding: utf-8 -*-
"""
Created on Tue May 26 02:47:08 2020

@author: dingxu
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture
import cv2
import os
import math
from itertools import combinations,permutations

fitsname1 = 'E:\\shunbianyuan\\newdata\\'+'20190603132646Auto.fit'
fitsname2 = 'E:\\shunbianyuan\\newdata\\'+'20190603132720Auto.fit'
#fitsname1 = 'E:\\shunbianyuan\\newdata\\'+'M31.fts'
#fitsname2 = 'E:\\shunbianyuan\\newdata\\'+'M31o.fts'
onehdu = fits.open(fitsname1)
imgdata1 = onehdu[0].data  #hdu[0].header
imgdata1 = np.rot90(imgdata1)

copydata1 = np.copy(imgdata1)
imgdata1 = np.float32(copydata1)
oneimgdata = imgdata1
timedate = onehdu[0].header['DATE']
#oneimgdata = signal.medfilt2d(imgdata1, kernel_size=5)  # 二维中值滤波
hang1,lie1 = oneimgdata.shape

twohdu = fits.open(fitsname2)
imgdata2 = twohdu[0].data  #hdu[0].header
#imgdata2 = np.rot90(imgdata2)
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
    daofind = DAOStarFinder(fwhm=12, threshold=5.*std)
    sources = daofind(img - median)

    #tezhen = np.transpose((sources['sharpness'], sources['roundness1'],sources['flux']))
    #print(sources[0])
    tezhen = np.transpose((sources['xcentroid'], sources['ycentroid']))
    posiandmag = np.transpose((sources['xcentroid'], sources['ycentroid'],sources['flux']))

    return tezhen,posiandmag.tolist()


###实现找星###
positions1,posiandmag1 =  findsource(oneimgdata)
positions2,posiandmag2 =  findsource(twoimgdata)

lenstar1 = len(positions1)
lenstar2 = len(positions2)


apertures1 = CircularAperture(positions1, r=13.)
apertures2 = CircularAperture(positions2, r=13.)

displayimage(oneimgdata,3,0)
apertures1.plot(color='blue', lw=1.5, alpha=0.5)

displayimage(twoimgdata,3,1)
apertures2.plot(color='blue', lw=1.5, alpha=0.5)

#newstar1 =  posiandmag1[np.lexsort(posiandmag1.T)]
#newstar12 =  posiandmag2[np.lexsort(posiandmag2.T)]

posiandmag1.sort(key=lambda x:x[2],reverse=True)
posiandmag2.sort(key=lambda x:x[2],reverse=True)

##选19颗亮星
lenstar = min(lenstar1,lenstar2)
lenstar = 14
posiandmag1 = posiandmag1[0:lenstar]
posiandmag2 = posiandmag2[0:lenstar]

sanjiao1 = list(combinations(posiandmag1,3))
sanjiao2 = list(combinations(posiandmag2,3))

def julisanjiao(sanjiao1,i):
    x1 = sanjiao1[i][0][0]
    y1 = sanjiao1[i][0][1]
    
    x2 = sanjiao1[i][1][0]
    y2 = sanjiao1[i][1][1]
    
    x3 = sanjiao1[i][2][0]
    y3 = sanjiao1[i][2][1]
    
    datadis1 = ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
    dS1S2 = math.sqrt(datadis1)
    
    datadis2 = ((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3))
    dS1S3 = math.sqrt(datadis2)
    
    datadis3 = ((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3))
    dS2S3 = math.sqrt(datadis3)
    
   #return [[x1,y1],[x2,y2],[x3,y3],[dS1S2,dS1S3,dS2S3]]
    return x1,y1,x2,y2,x3,y3,dS1S2,dS1S3,dS2S3


lensan1 = len(sanjiao1)
temp1 = np.zeros((lensan1,9),dtype = np.float32)
temp11 = np.zeros((lensan1,128),dtype = np.float32)
#keyimg2 = np.zeros((lenposition2,128),dtype = np.float32)

#lensan1 = len(sanjiao1)
#temp1 = []
for i in range (0,lensan1):
    jie1 = julisanjiao(sanjiao1,i)   
    if jie1[6]>jie1[7] and jie1[7]>jie1[8]:
        temp1[i,0:9] = jie1
        temp11[i,0:3] = temp1[i,6:9]

  
lensan2 = len(sanjiao2)    
temp2 = np.zeros((lensan1,9),dtype = np.float32)
temp22 = np.zeros((lensan1,128),dtype = np.float32)
for i in range (0,lensan2):
    jie2 = julisanjiao(sanjiao2,i) 
    if jie2[6]>jie2[7] and jie2[7]>jie2[8]:
        temp2[i,0:9] = jie2
        temp22[i,0:3] = temp2[i,6:9]

 
# FLANN 参数设计
FLANN_INDEX_KDTREE = 0
index_params = dict(algorithm = FLANN_INDEX_KDTREE, trees = 5)
search_params = dict(checks=50)
flann = cv2.FlannBasedMatcher(index_params,search_params)

matches = flann.knnMatch(temp11,temp22,k=2)

lenpipei = 0
knntemp1 = []
knntemp2 = []
for i, (m1, m2) in enumerate(matches):
    if m1.distance < 0.4 * m2.distance:# 两个特征向量之间的欧氏距离，越小表明匹配度越高。
        lenpipei = lenpipei+1
        knntemp1.append(m1.queryIdx)
        knntemp2.append(m1.trainIdx)

hmerge = np.hstack((imgdata1, imgdata2)) #水平拼接
displayimage(hmerge, 3, 2)

srckp1 = []
srckp2 = []
for i in range(lenpipei):
    x = knntemp1[i]
    y = knntemp2[i]
    x10 = temp1[x][0]
    y10 = temp1[x][1]
    x20 = temp1[x][2]
    y20 = temp1[x][3]
    x30 = temp1[x][4]
    y30 = temp1[x][5]
    
    srckp1.append(x10)
    srckp1.append(y10)
    srckp1.append(x20)
    srckp1.append(y20)
    srckp1.append(x30)
    srckp1.append(y30)
    src_pts = np.float32(srckp1).reshape(-1,2)
    
    x11 = temp2[y][0]
    y11 = temp2[y][1]
    x21 = temp2[y][2]
    y21 = temp2[y][3]
    x31 = temp2[y][4]
    y31 = temp2[y][5]
    
    srckp2.append(x11)
    srckp2.append(y11)
    srckp2.append(x21)
    srckp2.append(y21)
    srckp2.append(x31)
    srckp2.append(y31)
    dst_pts = np.float32(srckp2).reshape(-1,2)
    
    plt.plot([x10,x11+lie1],[y10,y11],linewidth = 0.8)
    plt.plot([x20,x21+lie1],[y20,y21],linewidth = 0.8)
    plt.plot([x30,x31+lie1],[y30,y31],linewidth = 0.8)
      

H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC,5.0)
newimg1 = cv2.warpPerspective(imgdata1, H, (lie1,hang1))

addimg = np.float32(newimg1) + np.float32(imgdata2)
minusimg = np.float32(newimg1) - np.float32(imgdata2)
displayimage(addimg, 3, 3)
displayimage(minusimg, 3, 4)
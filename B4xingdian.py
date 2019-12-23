# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 14:04:52 2019

@author: dingxu
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from skimage import measure

fitshdu = fits.open('E:\\BOOTES4\\20190606\\'+'201906061536370716.fit')
fitsdata = fitshdu[0].data

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

###求坐标##
def qiuzuobiao(img):

    threhold = np.mean(img)+3*np.std(img)    
    erzhiA1 = (img >= threhold)*1.0        
    labels = measure.label(erzhiA1,connectivity = 2) #8 连通区域标记
    print('regions number:',labels.max()+1) #显示连通区域块数(从 0 开始标

    regionnum = labels.max()+1

    hang,lie = img.shape
    plotx = np.zeros(regionnum,dtype = np.uint)
    ploty = np.zeros(regionnum,dtype = np.uint)

    for k in range(regionnum):
        sumx = 0
        sumy = 0
        area = 0
        for i in range(hang):
            for j in range(lie):
                if (labels[i][j] == k):
                    subimgthr = float(img[i,j])-float(threhold)
                    subimgthr0 = subimgthr if (subimgthr>0) else 0
                    sumx = sumx+i*subimgthr0
                    sumy = sumy+j*subimgthr0
                    area = area+subimgthr0
        try:
            plotx[k] = float(sumx/(area))
            ploty[k] = float(sumy/(area))
        except:
            plotx[k] = float(sumx/(area+0.0001))
            ploty[k] = float(sumy/(area+0.0001))
    return plotx,ploty,regionnum


mindata,maxdata = adjustimage(fitsdata, 3)
plotx,ploty,regionnum = qiuzuobiao(fitsdata)

plt.imshow(fitsdata, cmap='gray', vmin = mindata, vmax = maxdata)
for i in range(1,regionnum):
    plt.plot(ploty[i],plotx[i],'*')
plt.savefig('save.jpg')
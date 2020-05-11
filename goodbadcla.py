# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 18:23:33 2019

@author: dingxu
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
import xlwt

filename = 'J:\\bootes4\\xdr\\BOOTES-4\\2018\\20181118\\03088\\20181118132611-873-RA.fits'
fitshdu = fits.open(filename)
fitsdata = fitshdu[0].data

def guiyi(data):
    maxda = np.max(data)
    minda = np.min(data)
    guiyidata = (data-minda)/(maxda-minda)
    return guiyidata

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

    threhold = np.mean(img)+5*np.std(img)    
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
            plotx[k] = round(sumx/(area))
            ploty[k] = round(sumy/(area))
        except:
            plotx[k] = round(sumx/(area+0.0001))
            ploty[k] = round(sumy/(area+0.0001))
    return plotx,ploty,regionnum

def FWHMplot(x0,y0,width,img,i):
    x0 = int(x0)
    y0 = int(y0)
    pixlist = []
    for i in range((y0-width),(y0+width)):
        pixlist.append(img[x0,i])
    
    plt.figure(i)
    plt.plot(pixlist)
    plt.grid(color='r',linestyle='--',linewidth=2)

def calculateju(img,w,centerx,centery):
    count = w*w
    judata = np.copy(img[centerx-w:centerx+w,centery-w:centery+w])
    yijie = np.mean(judata)
    erjie = np.std(judata)
    newju = judata-yijie
    sanfangju = newju*newju*newju
    sumpingfang = np.sum(sanfangju)
    chun = sumpingfang/count
    sanjie = (chun)**(1/3)
    
    return yijie,erjie,sanjie

def writeexcel(hangcount,data_sheet,yijie,erjie,sanjie):
    data_sheet.write(hangcount, 0, yijie)
    data_sheet.write(hangcount, 1, erjie)
    data_sheet.write(hangcount, 2, sanjie)



guiyi = guiyi(fitsdata) 
mindata,maxdata = adjustimage(guiyi, 3)
#plotx,ploty,regionnum = qiuzuobiao(guiyi)


plt.figure(0)
plt.imshow(guiyi, cmap='gray', vmin = mindata, vmax = maxdata)

centery = 840
centerx = 634
hang = 0

plt.plot(centery,centerx,'*')
yijie,erjie,sanjie = calculateju(guiyi,20,centerx,centery)
FWHMplot(centerx,centery,20,guiyi,1)

sheet = 'B68'
workbook = xlwt.Workbook(encoding='utf-8')
data_sheet = workbook.add_sheet(sheet,cell_overwrite_ok = True)
writeexcel(hang,data_sheet,yijie,erjie,sanjie)
keepfilename =  'E:/B99.xls'
workbook.save(keepfilename)    
print(yijie,erjie,sanjie)


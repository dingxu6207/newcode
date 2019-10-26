# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 23:23:27 2019

@author: dingxu
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture
import math
import time
from photutils import create_matching_kernel
from photutils import TopHatWindow


start = time.time()

fitsname1 = 'd8549.0088.fits'
fitsname2 = 'd8549.0091.fits'
routename1 = 'E:\\shunbianyuan\\code\\phot\\'
routename2 = 'E:\\shunbianyuan\\code\\phot\\'

fitsname1 = routename1+fitsname1
fitsname2 = routename2+fitsname2


onehdu = fits.open(fitsname1)
oneimgdata = onehdu[1].data  #hdu[0].header

twohdu = fits.open(fitsname2)
twoimgdata = twohdu[1].data  #hdu[0].header

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
    daofind = DAOStarFinder(fwhm=4.0, threshold=5.*std)
    sources = daofind(img - median)

    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output
        #print(sources)

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    positionflux = np.transpose((sources['xcentroid'], sources['ycentroid'],  sources['flux']))
    mylist = positionflux.tolist()
    
    return sources,positions,mylist


###实现找星和算流量###
sources1,positions1,mylist1 =  findsource(oneimgdata)
mindata1,maxdata1 = adjustimage(oneimgdata,1)

sources2,positions2,mylist2 =  findsource(twoimgdata)
mindata2,maxdata2 = adjustimage(twoimgdata,1)
    
apertures1 = CircularAperture(positions1, r=10.)
apertures2 = CircularAperture(positions2, r=10.)

plt.figure(0)
plt.imshow(oneimgdata, cmap='gray', vmin = mindata1, vmax = maxdata1, origin='lower')
apertures1.plot(color='blue', lw=1.5, alpha=0.5)

plt.figure(1)
plt.imshow(twoimgdata, cmap='gray', vmin = mindata2, vmax = maxdata2, origin='lower')
apertures2.plot(color='blue', lw=1.5, alpha=0.5)


###匹配###
mylist1.sort(key=lambda x:x[2],reverse=True)
mylist2.sort(key=lambda x:x[2],reverse=True)

def suansanjiaoxing(listS1,listS2,listS3):
    duanchu = 0
    sumchen = 0
    x1 = float(listS1[1])
    y1 = float(listS1[0])
    x2 = float(listS2[1])
    y2 = float(listS2[0])
    x3 = float(listS3[1])
    y3 = float(listS3[0])
    
    datadis1 = ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
    dS1S2 = math.sqrt(datadis1)
    
    datadis2 = ((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3))
    dS1S3 = math.sqrt(datadis2)
    
    datadis3 = ((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3))
    dS2S3 = math.sqrt(datadis3)
        
    if (dS1S2 < dS1S3 and dS1S3 < dS2S3):
        duanchu = (dS2S3/dS1S2)
        sumchen = (x1-x3)*(x2-x3) + (y1-y3)*(y2-y3)
    
    if (dS1S2 < dS2S3 and dS2S3 < dS1S3):
        duanchu = (dS1S3/dS1S2)
        sumchen = (x1-x3)*(x2-x3) + (y1-y3)*(y2-y3)
        
    if (dS2S3 < dS1S3 and dS1S3 < dS1S2):
        duanchu = (dS1S2/dS2S3)
        sumchen = (x2-x1)*(x3-x1) + (y2-y1)*(y3-y1)    
        
    if (dS2S3 < dS1S2 and dS1S2 < dS1S3):
        duanchu = (dS1S3/dS2S3)
        sumchen = (x2-x1)*(x3-x1) + (y2-y1)*(y3-y1)    
           
    if (dS1S3 < dS1S2 and dS1S2 < dS2S3):
        duanchu = dS2S3/dS1S3
        sumchen = (x1-x2)*(x3-x2) + (y1-y2)*(y3-y2)
        
    if (dS1S3 < dS2S3 and dS2S3 < dS1S2):
        duanchu = dS1S2/dS1S3
        sumchen = (x1-x2)*(x3-x2) + (y1-y2)*(y3-y2)    
        
    return x1,y1,x2,y2,x3,y3,duanchu,sumchen

###第一张图算三角形###
lenlist1 = len(mylist1)
listsanjiao1 =  [0 for i in range(lenlist1)]
i = 6
for i in range(lenlist1):
    if (i <= (lenlist1-3)):
        l1x1,l1y1,l1x2,l1y2,l1x3,l1y3,l1duanchu,l1sumchen = suansanjiaoxing(mylist1[i],mylist1[i+1],mylist1[i+2])
        listsanjiao1[i] = (l1x1,l1y1,l1x2,l1y2,l1x3,l1y3,l1duanchu,l1sumchen)

###第二张图算三角形###
lenlist2 = len(mylist2)
listsanjiao2 =  [0 for i in range(lenlist2)]
i = 6
for i in range(lenlist2):
    if (i <= (lenlist2-3)):
        l2x1,l2y1,l2x2,l2y2,l2x3,l2y3,l2duanchu,l2sumchen = suansanjiaoxing(mylist2[i],mylist2[i+1],mylist2[i+2])
        listsanjiao2[i] = (l2x1,l2y1,l2x2,l2y2,l2x3,l2y3,l2duanchu,l2sumchen)
        

listtempA1 =  [0 for i in range(lenlist1)]
listtempA2 =  [0 for i in range(lenlist2)]
counttemp = 0
bilv = 0.1
dianchengsum = 200
#jiezhistar = min(A1regionnum, A2regionnum)
jiezhistar = 30     
for i in range(jiezhistar):
    for j in range(jiezhistar):
        if (abs(listsanjiao1[i][6]-listsanjiao2[j][6]) < bilv and abs(listsanjiao1[i][7]-listsanjiao2[j][7]) < dianchengsum):            
            listtempA1[counttemp] = listsanjiao1[i]
            listtempA2[counttemp] = listsanjiao2[j]
            counttemp = counttemp+1

###求平移###
def delhanshu(data0,data1,ydata0,ydata1):
    delx = ydata0-data0
    dely = ydata1-data1
        
    print('delx = ', delx)
    print('dely = ', dely)

    return delx,dely

data0 = listtempA1[0][0]  #x0
data1 = listtempA1[0][1]  #y0
data2 = listtempA1[0][2]  #x1
data3 = listtempA1[0][3]  #y1
data4 = listtempA1[0][4]  #x2
data5 = listtempA1[0][5]  #y2

###只比较x坐标###
###x0-x0###
abs01 = (listtempA1[0][0]-listtempA2[0][0]) - (listtempA1[0][2]-listtempA2[0][2])
abs02 = (listtempA1[0][0]-listtempA2[0][0]) - (listtempA1[0][2]-listtempA2[0][4])
###x0-x1###
abs10 = (listtempA1[0][0]-listtempA2[0][2]) - (listtempA1[0][2]-listtempA2[0][0])
abs12 = (listtempA1[0][0]-listtempA2[0][2]) - (listtempA1[0][2]-listtempA2[0][4])
###x0-x2###
abs20 = (listtempA1[0][0]-listtempA2[0][4]) - (listtempA1[0][2]-listtempA2[0][0])
abs21 = (listtempA1[0][0]-listtempA2[0][4]) - (listtempA1[0][2]-listtempA2[0][2])

redata = min(abs(abs01),abs(abs02),abs(abs10),abs(abs12),abs(abs20),abs(abs21))
if redata == abs(abs01):
    ydata0 = listtempA2[0][0] #x0
    ydata1 = listtempA2[0][1] #y0
    ydata2 = listtempA2[0][2] #x1
    ydata3 = listtempA2[0][3] #y1
    ydata4 = listtempA2[0][4] #x2
    ydata5 = listtempA2[0][5] #y2

if redata == abs(abs02):
    ydata0 = listtempA2[0][0] #x0
    ydata1 = listtempA2[0][1] #y0
    ydata2 = listtempA2[0][4] #x2
    ydata3 = listtempA2[0][5] #y2
    ydata4 = listtempA2[0][2] #x1
    ydata5 = listtempA2[0][3] #y1

if redata == abs(abs10):
    ydata0 = listtempA2[0][2] #x1
    ydata1 = listtempA2[0][3] #y1
    ydata2 = listtempA2[0][0] #x0
    ydata3 = listtempA2[0][1] #y0   
    ydata4 = listtempA2[0][4] #x2
    ydata5 = listtempA2[0][5] #y2

if redata == abs(abs12):
    ydata0 = listtempA2[0][2] #x1
    ydata1 = listtempA2[0][3] #y1
    ydata2 = listtempA2[0][4] #x2
    ydata3 = listtempA2[0][5] #y2
    ydata4 = listtempA2[0][0] #x0
    ydata5 = listtempA2[0][1] #y0
    

if redata == abs(abs20):
    ydata0 = listtempA2[0][4] #x2
    ydata1 = listtempA2[0][5] #y2
    ydata2 = listtempA2[0][0] #x0
    ydata3 = listtempA2[0][1] #y0
    ydata4 = listtempA2[0][2] #x1
    ydata5 = listtempA2[0][3] #y1

if redata == abs(abs21):
    ydata0 = listtempA2[0][4] #x2
    ydata1 = listtempA2[0][5] #y2
    ydata2 = listtempA2[0][2] #x1
    ydata3 = listtempA2[0][3] #y1 
    ydata4 = listtempA2[0][0] #x0
    ydata5 = listtempA2[0][1] #y0

delx1,dely1 =  delhanshu(data0,data1,ydata0,ydata1)   
delx2,dely2 =  delhanshu(data2,data3,ydata2,ydata3) 
delx3,dely3 =  delhanshu(data4,data5,ydata4,ydata5)
delx = round((delx1+delx2+delx3)/3)
dely = round((dely1+dely2+dely3)/3)

###画出匹配的点###
plt.figure(2)
plt.imshow(oneimgdata, cmap='gray', vmin = mindata1, vmax = maxdata1, origin='lower')
plt.plot(data1,data0,'*')
plt.plot(data3,data2,'*')
plt.plot(data5,data4,'*')

plt.figure(3)
plt.imshow(twoimgdata, cmap='gray', vmin = mindata2, vmax = maxdata2, origin='lower')
plt.plot(ydata1,ydata0,'*')
plt.plot(ydata3,ydata2,'*')
plt.plot(ydata5,ydata4,'*')

###图像平移###
###图像平移###
Ahang,Alie = oneimgdata.shape
newimage = np.zeros((Ahang,Alie),dtype = np.uint16)
if delx <= 0 and dely <= 0:
    newimage[0:Ahang+delx,0:Alie+dely] = oneimgdata[-delx:Ahang,-dely:Alie]
    
if delx <= 0 and dely >= 0:
    newimage[0:Ahang+delx,dely:Alie] = oneimgdata[-delx:Ahang,0:Alie-dely]
    
if delx >= 0 and dely >= 0:
    newimage[delx:Ahang,dely:Alie] = oneimgdata[0:Ahang-delx,0:Alie-dely]   
    
if delx >= 0 and dely <= 0:
    newimage[delx:Ahang,0:Alie+dely] = oneimgdata[0:Ahang-delx,-dely:Alie] 
    
jianimage = np.float32(newimage) - np.float32(twoimgdata)
jianimage = np.abs(jianimage)
minjian,maxjian = adjustimage(jianimage,1)
plt.figure(4)    
plt.imshow(jianimage,vmin=minjian,vmax=maxjian,cmap='gray') 
###PSF展宽###
def convimg(position,img1,img2,delx,dely):
    x = round(position[1])
    y = round(position[0])
    k = 15
    cropimg1 = img1[x-k:x+k,y-k:y+k]
    cropimg2 = img2[x-k-delx:x+k-delx,y-k-dely:y+k-dely]
    return cropimg1,cropimg2

cropimg1,cropimg2 = convimg(mylist2[38],twoimgdata,oneimgdata,delx,dely)
mincrop1,maxcrop1 = adjustimage(cropimg1,coffe = 1)
mincrop2,maxcrop2 = adjustimage(cropimg2,coffe = 1)
plt.figure(5)    
plt.imshow(cropimg1,vmin=mincrop1,vmax=maxcrop1,cmap='gray')

plt.figure(6)    
plt.imshow(cropimg2,vmin=mincrop2,vmax=maxcrop2,cmap='gray')

def make_blurred(gray, PSF, eps=0):
    input_fft = np.fft.fft2(gray)# 进行二维数组的傅里叶变换
    PSF_fft = np.fft.fft2(PSF)+ eps
    blurred = np.fft.ifft2(input_fft * PSF_fft)
    blurredimage = np.abs(np.fft.fftshift(blurred))
    
    return blurredimage

fcropimg1 = np.float64(cropimg1)
fcropimg2 = np.float64(cropimg2)
window = TopHatWindow(0.3)
kernel = create_matching_kernel(fcropimg1, fcropimg2, window)


bluimage = make_blurred(fcropimg1,kernel)

plt.figure(7)  
subimg = np.float64(bluimage) - np.float64(fcropimg2) 
subimg = np.abs(subimg)
plt.imshow(subimg,vmin=minjian,vmax=maxjian,cmap='gray')

plt.figure(8)  
fsubimg = np.float64(fcropimg1) - np.float64(fcropimg2) 
fsubimg = np.abs(fsubimg)
plt.imshow(fsubimg,vmin=minjian,vmax=maxjian,cmap='gray')

end = time.time()
print("运行时间:%.2f秒"%(end-start))

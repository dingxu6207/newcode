# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 00:22:23 2019

@author: dingxu
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture
import sanjiao
import time
from photutils import create_matching_kernel
from photutils import TopHatWindow
from photutils import CosineBellWindow
import copy


start = time.time()

fitsname1 = 'M31%20A6.fts'
fitsname2 = 'M31%20A6.fts'
routename1 = 'E:/shunbianyuan/ldf_download/20191030/'
routename2 = 'E:/shunbianyuan/ldf_download/20191031/'

fitsname1 = routename1+fitsname1
fitsname2 = routename2+fitsname2


onehdu = fits.open(fitsname1)
y1imgdata = onehdu[0].data  #hdu[0].header
oneimgdata = copy.deepcopy(y1imgdata)
oneimgdata = oneimgdata[300:600,300:600]

twohdu = fits.open(fitsname2)
y2imgdata = twohdu[0].data  #hdu[0].header
twoimgdata = copy.deepcopy(y2imgdata)
twoimgdata = twoimgdata[300:600,300:600]

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

listsanjiao1 = sanjiao.listsanjiao(mylist1)
listsanjiao2 = sanjiao.listsanjiao(mylist2)

listtempA1,listtempA2 = sanjiao.pipeisanjiao(listsanjiao1,listsanjiao2)

data = sanjiao.dianpipei(listtempA1[1],listtempA2[1])

###画出匹配的点###
plt.figure(2)
plt.imshow(oneimgdata, cmap='gray', vmin = mindata1, vmax = maxdata1, origin='lower')
plt.plot(data[1],data[0],'*')
plt.plot(data[3],data[2],'*')
plt.plot(data[5],data[4],'*')
#plt.plot(mylist2[8][0]+32,mylist2[8][1]+64,'*')


plt.figure(3)
plt.imshow(twoimgdata, cmap='gray', vmin = mindata2, vmax = maxdata2, origin='lower')
plt.plot(data[7],data[6],'*')
plt.plot(data[9],data[8],'*')
plt.plot(data[11],data[10],'*')


delx,dely = sanjiao.pingyiliang(data)
print(delx,dely)

###图像平移###
###图像平移###
Ahang,Alie = oneimgdata.shape
#newimage = np.zeros((Ahang,Alie),dtype = np.uint16)
if delx <= 0 and dely <= 0:
    newimage = np.zeros((Ahang+delx,Alie+dely),dtype = np.uint16)
    newimage2 = np.zeros((Ahang+delx,Alie+dely),dtype = np.uint16)
    newimage[0:Ahang+delx,0:Alie+dely] = oneimgdata[-delx:Ahang,-dely:Alie]
    newimage2[0:Ahang+delx,0:Alie+dely] = twoimgdata[0:Ahang+delx,0:Alie+dely]
    
if delx <= 0 and dely >= 0:
    newimage = np.zeros((Ahang+delx,Alie-dely),dtype = np.uint16)
    newimage2 = np.zeros((Ahang+delx,Alie+dely),dtype = np.uint16)
    newimage[0:Ahang+delx,dely:Alie] = oneimgdata[-delx:Ahang,0:Alie-dely]
    newimage2[0:Ahang+delx,0:Alie-dely] = twoimgdata[0:Ahang+delx,0:Alie-dely]
    
if delx >= 0 and dely >= 0:
    newimage = np.zeros((Ahang-delx,Alie-dely),dtype = np.uint16)
    newimage2 = np.zeros((Ahang-delx,Alie-dely),dtype = np.uint16)
    newimage = oneimgdata[0:Ahang-delx,0:Alie-dely]  
    newimage2 = twoimgdata[delx:Ahang,dely:Alie]
    
if delx >= 0 and dely <= 0:
    newimage = np.zeros((Ahang-delx,Alie+dely),dtype = np.uint16)
    newimage2 = np.zeros((Ahang-delx,Alie+dely),dtype = np.uint16)
    newimage[0:Ahang-delx,0:Alie+dely] = oneimgdata[0:Ahang-delx,-dely:Alie] 
    newimage2 = twoimgdata[0:Ahang-delx,0:Alie+dely]
    
jianimage = np.float32(newimage) - np.float32(newimage2)
#jianimage = np.abs(jianimage)
minjian,maxjian = adjustimage(jianimage,1)
plt.figure(4)    
plt.imshow(jianimage,vmin=minjian,vmax=maxjian,cmap='gray') 
###PSF展宽###
def make_blurred(gray, PSF):
    input_fft = np.fft.fft2(gray)# 进行二维数组的傅里叶变换
    PSF_fft = np.fft.fft2(PSF)
    blurred = np.fft.ifft2(input_fft * PSF_fft)
    blurredimage = np.abs(np.fft.fftshift(blurred))
    
    return blurredimage

newimage =np.float64(newimage)
newimage2 = np.float64(newimage2)
#window = TopHatWindow(0.3)
window = CosineBellWindow(alpha=0.35)
kernel = create_matching_kernel(newimage, newimage2,window)
bluimage = make_blurred(newimage,kernel)

plt.figure(5)  
subimg = np.float64(bluimage) - np.float64(newimage2) 
#subimg = np.abs(subimg)
minsub,maxsub = adjustimage(subimg,3)
plt.imshow(subimg,vmin=minsub,vmax=maxsub,cmap='gray')
filename = 'E:\\shunbianyuan\\BASS\\'+str(123)+'.jpg'
plt.savefig(filename)

plt.figure(6) 
minblu,maxblu = adjustimage(bluimage,1) 
plt.imshow(bluimage,vmin=minblu,vmax=maxblu,cmap='gray')
end = time.time()
print("运行时间:%.2f秒"%(end-start))
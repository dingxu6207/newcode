# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 15:57:14 2020

1、sessaligen
2、sessfindsatar
3、seseephometry
4、sessdisplay
@author: dingxu
"""

import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry


filetemp = []
count = 0
oripath = 'E:\\shunbianyuan\\dataxingtuan\\alngc7142\\'  #路径参数
for root, dirs, files in os.walk(oripath):
   for file in files:
       if (file[-5:] == '.fits'):
           count = count+1
           filetemp.append(file)
       
       
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
    #plt.plot(376.82, 137.232, '*')
    #plt.savefig(oripath+str(i)+'.jpg')

def photometryimg(positions, img, i):
    
    positionslist = positions.tolist()
    
    aperture = CircularAperture(positionslist, r=4) #2*FWHM
    annulus_aperture = CircularAnnulus(positionslist, r_in=8, r_out=12)#4-5*FWHM+2*FWHM
    apers = [aperture, annulus_aperture]
    
    displayimage(img, 1, i) ###画图1
    aperture.plot(color='blue', lw=0.5)
    annulus_aperture.plot(color='red', lw=0.2)
    plt.pause(0.001)
    plt.clf()
    
    phot_table = aperture_photometry(img, apers)
    bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
    bkg_sum = bkg_mean * aperture.area
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    phot_table['residual_aperture_sum'] = final_sum       
    posflux = np.column_stack((positions, phot_table['residual_aperture_sum']))  
    #return posflux
    magstar = 25 - 2.5*np.log10(abs(final_sum/1))
    return posflux,magstar

def sourcephotometry(targetx, targety, sumpho, threshold=10):
    hang,lie = sumpho.shape    
    for i in range(hang):
        delt = np.sqrt((targetx - sumpho[i][0])**2+(targety - sumpho[i][1])**2)
        if delt < threshold:
            #print(sumpho[i])
            mag = 25 - 2.5*np.log10(sumpho[i][2]/1) #90曝光时间
            #print(mag)
            return sumpho[i],mag

def pltquxian(datayuan):
    data = np.array(datayuan)  
    data1 = np.copy(data)
    u = np.mean(data1)   
    std = np.std(data1)
    error = data1[np.abs(data1 - u) > 3*std]
    data_c = data1[np.abs(data1 - u) <= 3*std] 
    print( len(error))
    return data_c
 
lacation = np.loadtxt(oripath+'location.txt')     

fitshdu = fits.open(oripath+filetemp[0])
fitsdata = fitshdu[0].data
displayimage(fitsdata, 1, 0)

jiaoyantemp = []
startemp = []

m = 2#行扫描 i = 39
n = 3 #列扫描 j = 39
for i in range(0, count):
    try:
        fitshdu = fits.open(oripath+filetemp[i])
        data = fitshdu[0].data    
        fitsdata = data[398*m:398+398*m,389*n:389+389*n]
        posflux,magstar = photometryimg(lacation, fitsdata, 1)           
        startemp.append(magstar) 
        arraytemp = np.array(startemp).T
        
        posflux1,mag1 = sourcephotometry(326, 134, posflux)  #比较星位置1
        posflux2,mag2 = sourcephotometry(361, 308, posflux)  #比较星位置2
              
        jiaoyan = mag1-mag2        
                
        jiaoyantemp.append(jiaoyan)
                
        print('ok')
    except:
        print('error!!!')
    
starlight = np.hstack((lacation, arraytemp)) 
np.savetxt('starlight.txt', starlight)   
#jiaoyandata = pltquxian(jiaoyan)       
plt.figure(2)
plt.plot(jiaoyantemp,'.')



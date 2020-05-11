# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 22:13:08 2019

@author: dingxu
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry

fitsname1 = 'E:\\BOOTES4\\20190606\\'+'201906061526400716.fit'

onehdu = fits.open(fitsname1)
imgdata1 = onehdu[0].data  #hdu[0].header
#obserdata = onehdu[0].header['DATE']

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
    plt.savefig(str(i)+'.jpg')

def findsource(img):    
    mean, median, std = sigma_clipped_stats(img, sigma=3.0) 
    daofind = DAOStarFinder(fwhm=10, threshold=5.*std)
    sources = daofind(img - median)
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    print(std)
    return positions

def photometryimg(positions, img, i):
    
    positionslist = positions.tolist()
    
    aperture = CircularAperture(positionslist, r=20) #2*FWHM
    annulus_aperture = CircularAnnulus(positionslist, r_in=50, r_out=70)#4-5*FWHM+2*FWHM
    apers = [aperture, annulus_aperture]
    
    displayimage(img,3,i) ###画图1
    aperture.plot(color='blue', lw=0.5)
    annulus_aperture.plot(color='red', lw=0.2)
    
    phot_table = aperture_photometry(img, apers)
    bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
    bkg_sum = bkg_mean * aperture.area
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    phot_table['residual_aperture_sum'] = final_sum       
    posflux = np.column_stack((positions, phot_table['residual_aperture_sum']))  
    return posflux

def sourcephotometry(targetx, targety, sumpho, threshold=10):
    hang,lie = sumpho.shape    
    for i in range(hang):
        delt = np.sqrt((targetx - sumpho[i][0])**2+(targety - sumpho[i][1])**2)
        if delt < threshold:
            print(sumpho[i])
            mag = 25 - 2.5*np.log10(sumpho[i][2]) #125717
            #print(mag)
            return sumpho[i],mag

def FWHMplot(x0,y0,width,imgdata1,i):
    x0 = int(x0)
    y0 = int(y0)
    pixlist = []
    for i in range((y0-width),(y0+width)):
        pixlist.append(imgdata1[x0,i])
    
    plt.figure(i)
    plt.plot(pixlist)
    plt.grid(color='r',linestyle='--',linewidth=2)
    
    
displayimage(imgdata1,3,0)  #画图0
positions1 = findsource(imgdata1)
apertures1 = CircularAperture(positions1, r=15.)
apertures1.plot(color='blue', lw=1.5, alpha=0.5)

posflux = photometryimg(positions1, imgdata1, 1)

posflux1,mag1 = sourcephotometry(570,230,posflux)
posflux2,mag2 = sourcephotometry(617,400,posflux)
posflux3,mag3 = sourcephotometry(748,411,posflux)
posflux4,mag4 = sourcephotometry(854,442,posflux)
posflux5,mag5 = sourcephotometry(160,336,posflux)
posflux6,mag6 = sourcephotometry(338,146,posflux)


print(mag1,mag2,mag3,mag4,mag5,mag6)
phtotemp = []
phtotemp.append(mag1)
phtotemp.append(mag2)
phtotemp.append(mag3)
phtotemp.append(mag4)
phtotemp.append(mag5)
phtotemp.append(mag6)

IRAF = [14.720,14.617,14.614,13.512,14.577,14.366] 

tuplematrix = np.where(posflux[:,2]==np.max(posflux[:,2]))
index = tuplematrix[0][0]
FWHMplot(posflux[index][1],posflux[index][0],15,imgdata1,2)
#FWHMplot(posflux2[1],posflux2[0],14,2) #画图2
        
plt.figure(3)
plt.plot(phtotemp,'--')
plt.plot(IRAF)

c = [phtotemp[i] - IRAF[i] for i in range(len(IRAF))]
plt.figure(4)
plt.plot(c)
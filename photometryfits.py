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

fitsname1 = 'E:\\BOOTES4\\03088\\'+'test.fits'

onehdu = fits.open(fitsname1)
imgdata1 = onehdu[0].data  #hdu[0].header

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
    daofind = DAOStarFinder(fwhm=5.5, threshold=4.*std)
    sources = daofind(img - median)
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

    return positions

def photometryimg(positions, img, i):
    
    positionslist = positions.tolist()
    
    aperture = CircularAperture(positionslist, r=11) #2*FWHM
    annulus_aperture = CircularAnnulus(positionslist, r_in=21, r_out=32)#4-5*FWHM+2*FWHM
    apers = [aperture, annulus_aperture]
    
    displayimage(img,3,i)
    aperture.plot(color='blue', lw=0.5)
    annulus_aperture.plot(color='red', lw=0.2)
    
    phot_table = aperture_photometry(img, apers)
    bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
    bkg_sum = bkg_mean * aperture.area
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    phot_table['residual_aperture_sum'] = final_sum       
    posflux = np.column_stack((positions, phot_table['residual_aperture_sum']))  
    return posflux

def sourcephotometry(targetx, targety, posflux, threshold=10):
    hang,lie = sumpho.shape    
    for i in range(hang):
        delt = np.sqrt((targetx - sumpho[i][0])**2+(targety - sumpho[i][1])**2)
        if delt < threshold:
            print(sumpho[i])
            mag = 25 - 2.5*np.log10(sumpho[i][2]) #125717
            print(mag)

displayimage(imgdata1,3,0)
positions1 = findsource(imgdata1)
apertures1 = CircularAperture(positions1, r=15.)
apertures1.plot(color='blue', lw=1.5, alpha=0.5)

sumpho = photometryimg(positions1, imgdata1, 1)

sourcephotometry(377.88,138.09,sumpho)



            
            
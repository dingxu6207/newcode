# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 18:57:49 2019

@author: dingxu
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture

filename = 'J:\\0716\\'+'201911292219550716.fit'
fitshdu = fits.open(filename)
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



def findsource(img):    
    mean, median, std = sigma_clipped_stats(img, sigma=3.0)
    daofind = DAOStarFinder(fwhm=8.5, threshold=5.*std)
    sources = daofind(img - median)

    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output
        #print(sources)

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    positionflux = np.transpose((sources['xcentroid'], sources['ycentroid'],  sources['flux']))
    #mylist = positionflux.tolist()
    
    return sources,positions,positionflux

def FWHMplot(x0,y0,width,img,i):
    x0 = int(x0)
    y0 = int(y0)
    pixlist = []
    for i in range((y0-width),(y0+width)):
        pixlist.append(img[x0,i])
    
    plt.figure(i)
    plt.plot(pixlist)
    plt.grid(color='r',linestyle='--',linewidth=2)
    
def displayimage(img, coff, i):
    minimg,maximg = adjustimage(img, coff)
    plt.figure(i)
    plt.imshow(img, cmap='gray', vmin = minimg, vmax = maximg)
    plt.savefig(str(i)+'.jpg')
    
sources1,positions1,mylist1 =  findsource(fitsdata)
apertures1 = CircularAperture(positions1, r=10.)
displayimage(fitsdata,3,0)
apertures1.plot(color='blue', lw=1.5, alpha=0.5)
tuplematrix = np.where(mylist1[:,2]==np.max(mylist1[:,2]))

index = tuplematrix[0][0]
FWHMplot(positions1[index][1],positions1[index][0],15,fitsdata,1)

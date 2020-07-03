# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 16:02:24 2020

@author: dingxu
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture
from scipy.optimize import curve_fit
from scipy import asarray as ar



#20190603132720Auto.fit
file  = '0.fits'
path = 'E:\\shunbianyuan\\todingx\\aligendata\\'
filename = path+file
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
    daofind = DAOStarFinder(fwhm = 5.06, threshold=5.*std)
    sources = daofind(img - median)

    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output
        #print(sources)

    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    positionflux = np.transpose((sources['xcentroid'], sources['ycentroid'],  sources['flux']))
    mylist = positionflux.tolist()
    
    return sources,positions,mylist

def FWHMplot(x0,y0,width,img,i):
    x0 = int(x0)
    y0 = int(y0)
    pixlist = []
    for i in range((y0-width),(y0+width)):
        pixlist.append(img[x0,i])
    
    plt.figure(i)
    plt.plot(pixlist)
    plt.grid(color='r',linestyle='--',linewidth=2)
    return pixlist

def gaussian(x,*param):
    return param[0]*np.exp(-np.power(x - param[1], 2.) / (2 * np.power(param[2], 2.)))+black

def displayimage(img, coff, i):
    minimg,maximg = adjustimage(img, coff)
    plt.figure(i)
    plt.imshow(img, cmap='gray', vmin = minimg, vmax = maximg)
    plt.savefig(str(i)+'.jpg')
    

sources1,positions1,mylist1 =  findsource(fitsdata)
apertures1 = CircularAperture(positions1, r=15.)
displayimage(fitsdata,1,0)
apertures1.plot(color='blue', lw=1.5, alpha=0.5)

np.savetxt(path+'location.txt', positions1,fmt='%f',delimiter=' ')

mylist1.sort(key=lambda x:x[2],reverse=True)
index = 3
width = 20

templist = FWHMplot(mylist1[index][1],mylist1[index][0],width,fitsdata,1)


x = ar(range(width*2))

black = np.mean(templist[0:10])

popt,pcov = curve_fit(gaussian,x,templist,p0=[3,4,3])
plt.plot(x,gaussian(x,*popt),'ro:',label='fit')
print(popt[1])
print('FWHM =',2.35482*popt[2])

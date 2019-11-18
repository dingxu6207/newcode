# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 14:33:59 2019

@author: dingxu
"""

import sim
from matplotlib import pyplot as plt
#import torch
#import torch.nn as nn
#from skimage.measure import compare_ssim
import glob
#import sim
import numpy as np
#from torch.utils import data
import random
import os

sav=[]
T0=sim.fitsread('E:/AST3/RA1626_DEC6300/L20170521_05700_162651+6300_60S_SI_00047.FITS')[0]
filelist=sorted(glob.glob('E:/AST3/RA1626_DEC6300/L2018*.FITS'))
#T0=(T0[1000:8000,1000:8000])
#T0=np.sqrt(T0)
T0-=np.median(T0)
sav.append(T0)
for file in filelist:
    T=sim.fitsread(file)[0]
#    T=np.sqrt(T)
#    T=(T[1000:8000,1000:8000])
    T-=np.median(T)
    sav.append(T)
sav=np.array(sav)
#sim.fitswrite('Z1.fits',sav,header=None)    
#
#Z=sim.fitsread('Z1.fits')[0]
Z=np.array(sav)
dx0,dy0=0,0
for i in range(1,4):
    dx,dy,cor=sim.xcorrcenter(Z[0][4000:6000,4000:6000],Z[i][4000:6000,4000:6000])
    print(i,dx,dy)
    Z[i]=sim.immove(Z[i],dx,dy)
    dx0=np.max((dx,dx0))
    dy0=np.max((dy,dy0))
    print(i)
Z0=Z[0]
Z1=np.mean(Z[1:3,:,:],axis=0)    
#Z1=Z[2]
dx0=abs(np.round(dx0).astype(int))
dy0=abs(np.round(dy0).astype(int))
    
Z0=Z0[dx0:-dx0,dy0:-dy0]
Z1=Z1[dx0:-dx0,dy0:-dy0]
Z=np.stack((Z0,Z1))
plt.figure()
sim.showim(Z1/Z0)
plt.figure()
sim.showim(Z1)

#sim.fitswrite('data.fits',Z)
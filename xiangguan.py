# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 23:10:29 2020

@author: dingxu
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fft
import cv2 

fitsname1 = 'E:\\shunbianyuan\\dataxingtuan\\ngc7142\\d4738777L016m001.fit'
#fitsname1 = 'E:\\shunbianyuan\\todingx\\origindata\\ftboYFAk120299.fits'
onehdu = fits.open(fitsname1)
imgdata1 = onehdu[0].data  #hdu[0].header
copydata1 = np.copy(imgdata1)
onedata = np.float32(copydata1)

fitsname2 = 'E:\\shunbianyuan\\dataxingtuan\\ngc7142\\d4738777L016m008.fit'
#fitsname2 = 'E:\\shunbianyuan\\todingx\\origindata\\ftboYFBb250086.fits'
onehdu2 = fits.open(fitsname2)
imgdata2 = onehdu2[0].data  #hdu[0].header
copydata2 = np.copy(imgdata2)
twodata = np.float32(copydata2)

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
    plt.savefig(str(i)+'.jpg')

def zscore2(im):
    im = (im - np.mean(im)) / im.std()
    return im

def xcorrcenter(standimage, compimage, R0=3, flag=0):
    # if flag==1,standimage 是FFT以后的图像，这是为了简化整数象元迭代的运算量。直接输入FFT以后的结果，不用每次都重复计算
    try:
        M, N = standimage.shape

        standimage = zscore2(standimage)
        s = fft.fft2(standimage)

        compimage = zscore2(compimage)
        c = np.fft.ifft2(compimage)

        sc = s * c
        im = np.abs(fft.fftshift(fft.ifft2(sc)))  # /(M*N-1);%./(1+w1.^2);
        cor = im.max()
        if cor == 0:
            return 0, 0, 0

        M0, N0 = np.where(im == cor)
        m, n = M0[0], N0[0]

        if flag:
            m -= M / 2
            n -= N / 2
            # 判断图像尺寸的奇偶
            if np.mod(M, 2): m += 0.5
            if np.mod(N, 2): n += 0.5

            return m, n, cor
        # 求顶点周围区域的最小值
        immin = im[(m - R0):(m + R0 + 1), (n - R0):(n + R0 + 1)].min()
        # 减去最小值
        im = np.maximum(im - immin, 0)
        # 计算重心
        x, y = np.mgrid[:M, :N]
        area = im.sum()
        m = (np.double(im) * x).sum() / area
        n = (np.double(im) * y).sum() / area
        # 归算到原始图像
        m -= M / 2
        n -= N / 2
        # 判断图像尺寸的奇偶
        if np.mod(M, 2): m += 0.5
        if np.mod(N, 2): n += 0.5
    except:
        print('Err in align_Subpix routine!')
        m, n, cor = 0, 0, 0
    return m, n, cor



def immove(image, dx, dy):
    """
    image shift by subpix
    """
    # The shift corresponds to the pixel offset relative to the reference image
    from scipy.ndimage import fourier_shift
    if dx == 0 and dy == 0:
        offset_image = image
    else:
        shift = (dx, dy)
        offset_image = fourier_shift(np.fft.fft2(image), shift)
        offset_image = np.real(fft.ifft2(offset_image))

    return offset_image

rows, cols = copydata1.shape
dx,dy,cor=xcorrcenter(copydata1, copydata2)
#print(dx, dy)
#print(dx,dy)
M = np.float32([[1,0,dx],[0,1,dy]])
newdata = cv2.warpAffine(twodata, M, (cols,rows))
#newdata = immove(imgdata2,dx,dy)
displayimage(copydata1,1,0)
displayimage(copydata2,1,1)
displayimage(newdata,1,2)
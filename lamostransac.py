# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:28:15 2020

@author: dingxu
"""

import os
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
from scipy import signal
import operator as op

zuoxian = 2569-60
youxian = 2597+60

path = 'E:/pytest/B6202/'
strfile = 'E:/pytest/B6201/spec-55862-B6202_sp01-026.fits'
phdulist = fits.open(strfile)

#print(phdulist.info)
#print(phdulist[0].header)

flux = phdulist[0].data[0]
Liwave = phdulist[0].data[2]

duanwave = Liwave[zuoxian:youxian]
duanflux = flux[zuoxian:youxian]

n_dot = youxian-zuoxian
n_order = 3

x = np.linspace(0,1,n_dot) 

p = np.poly1d(np.polyfit(x,duanflux,n_order)) 

plt.figure(0)
plt.plot(Liwave,flux)


plt.figure(1)
plt.plot(duanwave,duanflux,duanwave,p(x))

plt.figure(2)
m = duanflux/p(x)
plt.plot(duanwave,m)


n_order = 3

sample = np.vstack((x,duanflux)).T
plt.figure(3)
plt.plot(x,duanflux)

class Ransac:
    a = 0.
    b = 0.
    c = 0.
    d = 0.
    n_order = 3  #ax**3+bx**2+cx+d
    def least_square(self,samples):
        ##最小二乘法
        x = samples[:,0]
        #x = np.linspace(0, 1, n_dot)
        y = samples[:,1]
        p = np.poly1d(np.polyfit(x,y,n_order)) 
        
        #print(p.coeffs)
        a,b,c,d = p.coeffs
        
        return a,b,c,d

    def isRepeat(self,sour,tar):
        #判断是否含有重复样本
        for i in range(len(sour)):
            if (op.eq(list(sour[i]), list(tar))):
                    return True
        return False

    def random_samples(self,samples,points_ratio):
        ## 随机采样（无重复样本）
        number = len(samples)
        inliers_num = int(number * points_ratio)
        inliers = []
        outliers = []
        cur_num = 0
        while cur_num != inliers_num:
            seed = np.random.randint(0,number)
            sap_cur = samples[seed]
            if not self.isRepeat(inliers,sap_cur):
                cur_num = cur_num +1
                inliers.append(list(sap_cur))
        for i in range(number):
            if not self.isRepeat(inliers,samples[i]):
                outliers.append(list(samples[i]))
        return np.array(inliers),np.array(outliers)

    def fun_plot(self,sample,a,b,c,d):
        data_x = np.linspace(0, 1, n_dot)
        data_y = [a * x**3 + b*x**2+c*x+d for x in data_x]
        plt.figure(4)
        plt.ion()
        plt.plot(data_x,data_y,'r')
        plt.plot(sample[:,0],sample[:,1])
        plt.show()
        plt.pause(0.05)
        plt.clf()

    def ransac(self,samples, points_ratio = 0.4, epoch = 100, reject_dis = 2 ,inliers_ratio = 0.7):
        # samples 输入样本，形如 [[x1 ,yi],[x2, y2]]
        # point_ratio  随机选择样本点的比例
        # epoch    迭代轮数
        # reject_dis  小于此阈值将outliers加入inliers
        # inliers_ratio  有效inliers最低比例

        inliers_num_cur = 0
        for i in range(epoch):
            inliers,outliers = self.random_samples(samples,points_ratio)
            #print(inliers.shape)
            a,b,c,d = self.least_square(inliers)
            # self.fun_plot(samples,weight_cur,bias_cur)
            for j in range(len(outliers)):
                distance = np.abs(a*(outliers[j,0]**3)+b*(outliers[j,0]**2)+c*(outliers[j,0])+d - outliers[j,1])
                if distance <=  reject_dis:
                    inliers = np.vstack((inliers,outliers[j]))
                    
            a,b,c,d = self.least_square(inliers)
            self.fun_plot(samples,a,b,c,d)
            
            if len(inliers) >= len(samples)* inliers_ratio:
               if len(inliers) > inliers_num_cur:
                    self.a = a
                    self.b = b
                    self.c = c
                    self.d = d                   
                    inliers_num_cur = len(inliers)
                    print(i,len(inliers))
        
test = Ransac()
test.ransac(sample)
data_x = np.linspace(0,1,n_dot)
data_y = [test.a * x**3 + test.b*x**2+test.c*x+test.d for x in data_x]
plt.plot(duanwave, sample[:, 1])
plt.plot(duanwave,data_y,'r')
plt.show()
plt.pause(4)

plt.figure(5)
guiyi = sample[:, 1]/data_y
plt.plot(duanwave, guiyi)


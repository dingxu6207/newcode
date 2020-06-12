# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 14:07:10 2020

@author: dingxu
"""

import os
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
from scipy import signal
import operator as op


n_order = 4

matrix = np.loadtxt('pu1.txt')

duanwave  = matrix[:,0]
duanflux = matrix[:,1]

duanwave = duanwave[258:1427]
duanflux = duanflux[258:1427]

x = (duanwave-np.min(duanwave))/(np.max(duanwave)-np.min(duanwave))
p = np.poly1d(np.polyfit(x,duanflux,n_order))

plt.figure(0)
plt.plot(duanwave,duanflux)
plt.plot(duanwave,p(x))

sample = np.vstack((x,duanflux)).T
plt.figure(1)
plt.plot(x,duanflux)

class Ransac:
    a = 0.
    b = 0.
    c = 0.
    d = 0.
    e = 0
    n_order = 4  #ax**3+bx**2+cx+d
    def least_square(self,samples):
        ##最小二乘法
        x = samples[:,0]
        #x = np.linspace(0, 1, n_dot)
        y = samples[:,1]
        p = np.poly1d(np.polyfit(x,y,n_order)) 
        
        #print(p.coeffs)
        a,b,c,d,e = p.coeffs
        
        return a,b,c,d,e

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

    def fun_plot(self,sample,a,b,c,d,e):
        #data_x = np.linspace(0, 1, n_dot)
        data_x = sample[:,0]
        data_y = [a * x**4 + b*x**3+c*x**2+d*x+e for x in data_x]
        plt.figure(2)
        plt.ion()
        plt.plot(data_x,data_y,'r')
        plt.plot(sample[:,0],sample[:,1])
        plt.show()
        plt.pause(0.05)
        plt.clf()

    def ransac(self,samples, points_ratio = 0.005, epoch = 1000, reject_dis = 0.02e-15 ,inliers_ratio = 0.1):
        # samples 输入样本，形如 [[x1 ,yi],[x2, y2]]
        # point_ratio  随机选择样本点的比例
        # epoch    迭代轮数
        # reject_dis  小于此阈值将outliers加入inliers
        # inliers_ratio  有效inliers最低比例

        inliers_num_cur = 0
        for i in range(epoch):
            inliers,outliers = self.random_samples(samples,points_ratio)
            #print(inliers.shape)
            a,b,c,d,e = self.least_square(inliers)
            # self.fun_plot(samples,weight_cur,bias_cur)
            for j in range(len(outliers)):
                distance = np.abs(a*(outliers[j,0]**4)+b*(outliers[j,0]**3)+c*(outliers[j,0]**2)+d*outliers[j,0]+e - outliers[j,1])
                if distance <=  reject_dis:
                    inliers = np.vstack((inliers,outliers[j]))
                    
            a,b,c,d,e = self.least_square(inliers)
            self.fun_plot(samples,a,b,c,d,e)
                        
            if len(inliers) >= len(samples)* inliers_ratio:
               if len(inliers) > inliers_num_cur:
                    self.a = a
                    self.b = b
                    self.c = c
                    self.d = d 
                    self.e = e   
                    self.inliers = inliers
                    inliers_num_cur = len(inliers)
                    print(i,len(inliers))
        
test = Ransac()
test.ransac(sample)
data_x = sample[:,0]
data_y = [test.a * x**4 + test.b*x**3+test.c*x**2+test.d*x+test.e for x in data_x]
plt.plot(duanwave, sample[:, 1])
plt.plot(duanwave,data_y,'r')

induanwave = test.inliers[:,0]*((np.max(duanwave)-np.min(duanwave)))+np.min(duanwave)
plt.plot(induanwave,test.inliers[:,1],'.')
plt.show()
plt.pause(4)

plt.figure(3)
guiyi = sample[:, 1]/data_y
plt.plot(duanwave, guiyi)
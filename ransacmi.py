# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 21:07:42 2020

@author: dingxu
"""

import os
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
from scipy import signal
import operator as op
from scipy.optimize import curve_fit

n_order = 3

matrix = np.loadtxt('pu1.txt')

duanwave  = matrix[:,0]
duanflux = matrix[:,1]

duanwave = duanwave[220:1427]
duanflux = duanflux[220:1427]

def fund(x, a, b):
    return a*(x**b )

x = duanwave/5100
popt, pcov = curve_fit(fund, x, duanflux)

ydata = popt[0]*x**popt[1]
plt.figure(0)
plt.plot(duanwave,ydata)
plt.plot(duanwave,duanflux)
sample = np.vstack((x,duanflux)).T

class Ransac:
    a = 0.
    b = 0.
    
    def fund(self,x, a, b):
        return a*(x**b )
    
    
    def least_square(self,samples):
        x = samples[:,0]
        duanflux = samples[:,1]
        popt, pcov = curve_fit(self.fund, x, duanflux)
        a = popt[0]
        b = popt[1]
        return a,b

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

    def fun_plot(self,sample,a,b):
        data_x = sample[:,0] 
        data_y = a*(data_x**b)
        plt.figure(1)
        plt.ion()
        plt.plot(data_x,data_y,'r')
        plt.plot(sample[:,0],sample[:,1])
        plt.show()
        plt.pause(0.05)
        plt.clf()

    def ransac(self,samples, points_ratio = 0.003, epoch = 400, reject_dis = 0.02e-15 ,inliers_ratio = 0.1):
        # samples 输入样本，形如 [[x1 ,yi],[x2, y2]]
        # point_ratio  随机选择样本点的比例
        # epoch    迭代轮数
        # reject_dis  小于此阈值将outliers加入inliers
        # inliers_ratio  有效inliers最低比例

        inliers_num_cur = 0
        for i in range(epoch):
            inliers,outliers = self.random_samples(samples,points_ratio)
            #print(inliers.shape)
            a,b = self.least_square(inliers)
            # self.fun_plot(samples,weight_cur,bias_cur)
            for j in range(len(outliers)):
                distance = np.abs(a*(outliers[j,0]**b) - outliers[j,1])
                if distance <=  reject_dis:
                    inliers = np.vstack((inliers,outliers[j]))
                    
            a,b = self.least_square(inliers)
            self.fun_plot(samples,a,b)
            
            if len(inliers) >= len(samples)* inliers_ratio:
               if len(inliers) > inliers_num_cur:
                    self.a = a
                    self.b = b     
                    self.inliers = inliers
                    
                    inliers_num_cur = len(inliers)
                    print(i,len(inliers))

test = Ransac()
test.ransac(sample)
data_x = sample[:,0]
data_y = [test.a * x**test.b for x in data_x]
plt.plot(duanwave, sample[:, 1])
plt.plot(duanwave,data_y,'r')
plt.show()
plt.pause(4)

induanwave = test.inliers[:,0]*5100
plt.plot(induanwave,test.inliers[:,1],'.')
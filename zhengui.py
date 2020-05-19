# -*- coding: utf-8 -*-
"""
Created on Tue May 19 15:38:49 2020

@author: dingxu
"""


import numpy as np
import cv2

a=np.mat([[1,1,1],[1,3,-1],[2,5,2],[3,-1,5],])#系数矩阵
b=np.mat([2,-1,1,-2]).T    #常数项列矩阵


m = np.dot(a.T,a)
n = np.dot(a.T,b)
mli = np.linalg.inv(m)
x = np.dot(mli,n)

positions1 = np.float32([[713.735,540.193],[154.701,297.219],[334.049,107.978]])

positions2 = np.float32([[712.337,538.53],[153.574,295.828],[332.549,106.124]])


post1 = np.float32([[713.735,540.193]])
post2 = np.float32([[712.337,538.53]])

M = cv2.getAffineTransform(positions1,positions2)

new1 = post1[0][0]*M[0][0]+post1[0][1]*M[0][1]+M[0][2]
new2 = post1[0][0]*M[1][0]+post1[0][1]*M[1][1]+M[1][2]
print(new1)
print(new2)

juli = np.sqrt((new1-post2[0][0])**2+(new2-post2[0][1])**2)

print(juli)

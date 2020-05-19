# -*- coding: utf-8 -*-
"""
Created on Tue May 19 22:08:44 2020

@author: dingxu
"""

import numpy as np
import matplotlib.pyplot as plt

#a=np.mat([[1,1,1],[1,3,-1],[2,5,2],[3,-1,5],])#系数矩阵
a = np.float32([[334.049,107.978,1,0,0,0],
            [0,0,0,334.049,107.978,1],
            
            [960.041,113.597,1,0,0,0],
            [0,0,0,960.041,113.597,1],
            
            [567.15,187.546,1,0,0,0],
            [0,0,0,567.15,187.546,1],
            
            [154.701,297.219,1,0,0,0],
            [0,0,0,154.701,297.219,1],
            
            [609.576,355.069,1,0,0,0],
            [0,0,0,609.576,355.069,1],
            
            [739.827,366.559,1,0,0,0],
            [0,0,0,739.827,366.559,1],
            
            [845.845,399.465,1,0,0,0],
            [0,0,0,845.845,399.465,1],
            
            [713.735,540.193,1,0,0,0],
            [0,0,0,713.735,540.193,1],
            
            [643.34,659.523,1,0,0,0],
            [0,0,0,643.34,659.523,1],
            
            [632.83,673.916,1,0,0,0],
            [0,0,0,632.83,673.916,1],
            
            [914.795,734.975,1,0,0,0],
            [0,0,0,914.795,734.975,1],
            
            [491.198,940.262,1,0,0,0],
            [0,0,0,491.198,940.262,1],
            
            ])

b=np.float32([332.549,
          106.124,
          
          958.312,
          112.69,
          
          566.238,
          186.695,
          
          153.574,
          295.828,
          
          607.67,
          353.778,
          
          738.635,
          365.537,
          
          844.78,
          398.529,
          
          712.337,
          538.53,
          
          641.752,
          657.18,
          
          631.345,
          671.994,
          
          913.196,
          732.829,
          
          489.609,
          938.567
          
          ]).T    #常数项列矩阵


m = np.dot(a.T,a)
n = np.dot(a.T,b)
mli = np.linalg.inv(m)
x = np.dot(mli,n)


print(a[0][0])
newx = a[0][0]*x[0]+a[0][1]*x[1]+x[2]
newy = a[1][3]*x[3]+a[1][4]*x[4]+x[5]

print(newx,newy)

juli = np.sqrt((newx-b[0])**2+(newy-b[1])**2)

temp = []
for i in range (12):
    p = 2*i
    pstx = a[p][0]*x[0]+a[p][1]*x[1]+x[2]
    jux = pstx-b[p]
    
    p = p+1
    psty = a[p][3]*x[3]+a[p][4]*x[4]+x[5]
    juy = psty -b[p]
    
    juli = np.sqrt(jux**2+juy**2)
    temp.append(juli)
    
plt.plot(temp)
plt.title('6 constant model')  

print(np.mean(temp))
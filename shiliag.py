# -*- coding: utf-8 -*-
"""
Created on Mon May 18 21:31:22 2020

@author: dingxu
"""
import matplotlib.pyplot as plt

tx0 = 334.049
ty0 = 107.978

tx1 = 332.549
ty1 = 106.124

plt.xlim(0, 400)
plt.ylim(0, 400)
plt.annotate(' ',xy=(tx0,ty0),xytext=(tx1,ty1),arrowprops=dict(arrowstyle="->",connectionstyle="arc3"))

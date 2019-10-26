# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 18:50:30 2019

@author: dingxu
"""

import os

i = 0
for i in range(60) :
    cmd = 'python phot.py'+' '+str(i)
    try:
        os.system(cmd)
    except:
        continue
        
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 00:17:16 2019

@author: dingxu
"""
import os
listroute = [0 for i in range(8)]
listroute[0] = 'http://psp.china-vo.org/pspdata/2019/11/20191104/M31-18object/'
listroute[1] = 'http://psp.china-vo.org/pspdata/2019/11/20191102/M31-18object/'
listroute[2] = 'http://psp.china-vo.org/pspdata/2019/11/20191101/M31-18object/'
listroute[3] = 'http://psp.china-vo.org/pspdata/2019/10/20191031/M31-18object/'
listroute[4] = 'http://psp.china-vo.org/pspdata/2019/10/20191030/M31-18object/'
listroute[5] = 'http://psp.china-vo.org/pspdata/2019/10/20191029/M31-18object/'
listroute[6] = 'http://psp.china-vo.org/pspdata/2019/10/20191028/M31-18object/'
listroute[7] = 'http://psp.china-vo.org/pspdata/2019/10/20191027/M31-18object/'

listpath = [0 for i in range(8)]
listpath[0] = 'E:/shunbianyuan/ldf_download/20191104'
listpath[1] = 'E:/shunbianyuan/ldf_download/20191102'
listpath[2] = 'E:/shunbianyuan/ldf_download/20191101'
listpath[3] = 'E:/shunbianyuan/ldf_download/20191031'
listpath[4] = 'E:/shunbianyuan/ldf_download/20191030'
listpath[5] = 'E:/shunbianyuan/ldf_download/20191029'
listpath[6] = 'E:/shunbianyuan/ldf_download/20191028'
listpath[7] = 'E:/shunbianyuan/ldf_download/20191027'

for i in range(8):
    cmd = 'python downdata.py'+' '+listroute[i]+' '+listpath[i]
    os.system(cmd)
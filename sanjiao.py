# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 00:00:45 2019

@author: dingxu
"""
import math

def suansanjiaoxing(listS1,listS2,listS3):
    duanchu = 0
    sumchen = 0
    x1 = float(listS1[1])
    y1 = float(listS1[0])
    x2 = float(listS2[1])
    y2 = float(listS2[0])
    x3 = float(listS3[1])
    y3 = float(listS3[0])
    
    datadis1 = ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
    dS1S2 = math.sqrt(datadis1)
    
    datadis2 = ((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3))
    dS1S3 = math.sqrt(datadis2)
    
    datadis3 = ((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3))
    dS2S3 = math.sqrt(datadis3)
        
    if (dS1S2 < dS1S3 and dS1S3 < dS2S3):
        duanchu = (dS2S3/dS1S2)
        sumchen = (x1-x3)*(x2-x3) + (y1-y3)*(y2-y3)
    
    if (dS1S2 < dS2S3 and dS2S3 < dS1S3):
        duanchu = (dS1S3/dS1S2)
        sumchen = (x1-x3)*(x2-x3) + (y1-y3)*(y2-y3)
        
    if (dS2S3 < dS1S3 and dS1S3 < dS1S2):
        duanchu = (dS1S2/dS2S3)
        sumchen = (x2-x1)*(x3-x1) + (y2-y1)*(y3-y1)    
        
    if (dS2S3 < dS1S2 and dS1S2 < dS1S3):
        duanchu = (dS1S3/dS2S3)
        sumchen = (x2-x1)*(x3-x1) + (y2-y1)*(y3-y1)    
           
    if (dS1S3 < dS1S2 and dS1S2 < dS2S3):
        duanchu = dS2S3/dS1S3
        sumchen = (x1-x2)*(x3-x2) + (y1-y2)*(y3-y2)
        
    if (dS1S3 < dS2S3 and dS2S3 < dS1S2):
        duanchu = dS1S2/dS1S3
        sumchen = (x1-x2)*(x3-x2) + (y1-y2)*(y3-y2)    
        
    return x1,y1,x2,y2,x3,y3,duanchu,sumchen

def listsanjiao(mylist1):
    mylist1.sort(key=lambda x:x[2],reverse=True)
    lenlist1 = len(mylist1)
    listsanjiao1 =  [0 for i in range(lenlist1)]
    i = 6
    for i in range(lenlist1):
        if (i <= (lenlist1-3)):
            l1x1,l1y1,l1x2,l1y2,l1x3,l1y3,l1duanchu,l1sumchen = suansanjiaoxing(mylist1[i],mylist1[i+1],mylist1[i+2])
            listsanjiao1[i] = (l1x1,l1y1,l1x2,l1y2,l1x3,l1y3,l1duanchu,l1sumchen)
    return listsanjiao1

def pipeisanjiao(listsanjiao1,listsanjiao2,bilv = 0.1,dianchengsum = 200,jiezhistar = 30 ):
    lenlist1 = len(listsanjiao1)
    lenlist2 = len(listsanjiao2)
    listtempA1 =  [0 for i in range(lenlist1)]
    listtempA2 =  [0 for i in range(lenlist2)]
    counttemp = 0
    for i in range(jiezhistar):
        for j in range(jiezhistar):
            if (abs(listsanjiao1[i][6]-listsanjiao2[j][6]) < bilv and abs(listsanjiao1[i][7]-listsanjiao2[j][7]) < dianchengsum):            
                listtempA1[counttemp] = listsanjiao1[i]
                listtempA2[counttemp] = listsanjiao2[j]
                counttemp = counttemp+1
    return listtempA1,listtempA2

def dianpipei(listtempA1,listtempA2):
    data0 = listtempA1[0]  #x0
    data1 = listtempA1[1]  #y0
    data2 = listtempA1[2]  #x1
    data3 = listtempA1[3]  #y1
    data4 = listtempA1[4]  #x2
    data5 = listtempA1[5]  #y2

###只比较x坐标###
###x0-x0###
    abs01 = (listtempA1[0]-listtempA2[0]) - (listtempA1[2]-listtempA2[2])
    abs02 = (listtempA1[0]-listtempA2[0]) - (listtempA1[2]-listtempA2[4])
###x0-x1###
    abs10 = (listtempA1[0]-listtempA2[2]) - (listtempA1[2]-listtempA2[0])
    abs12 = (listtempA1[0]-listtempA2[2]) - (listtempA1[2]-listtempA2[4])
    ###x0-x2###
    abs20 = (listtempA1[0]-listtempA2[4]) - (listtempA1[2]-listtempA2[0])
    abs21 = (listtempA1[0]-listtempA2[4]) - (listtempA1[2]-listtempA2[2])

    redata = min(abs(abs01),abs(abs02),abs(abs10),abs(abs12),abs(abs20),abs(abs21))
    if redata == abs(abs01):
        ydata0 = listtempA2[0] #x0
        ydata1 = listtempA2[1] #y0
        ydata2 = listtempA2[2] #x1
        ydata3 = listtempA2[3] #y1
        ydata4 = listtempA2[4] #x2
        ydata5 = listtempA2[5] #y2

    if redata == abs(abs02):
        ydata0 = listtempA2[0] #x0
        ydata1 = listtempA2[1] #y0
        ydata2 = listtempA2[4] #x2
        ydata3 = listtempA2[5] #y2
        ydata4 = listtempA2[2] #x1
        ydata5 = listtempA2[3] #y1

    if redata == abs(abs10):
        ydata0 = listtempA2[2] #x1
        ydata1 = listtempA2[3] #y1
        ydata2 = listtempA2[0] #x0
        ydata3 = listtempA2[1] #y0   
        ydata4 = listtempA2[4] #x2
        ydata5 = listtempA2[5] #y2

    if redata == abs(abs12):
        ydata0 = listtempA2[2] #x1
        ydata1 = listtempA2[3] #y1
        ydata2 = listtempA2[4] #x2
        ydata3 = listtempA2[5] #y2
        ydata4 = listtempA2[0] #x0
        ydata5 = listtempA2[1] #y0
    

    if redata == abs(abs20):
        ydata0 = listtempA2[4] #x2
        ydata1 = listtempA2[5] #y2
        ydata2 = listtempA2[0] #x0
        ydata3 = listtempA2[1] #y0
        ydata4 = listtempA2[2] #x1
        ydata5 = listtempA2[3] #y1

    if redata == abs(abs21):
        ydata0 = listtempA2[4] #x2
        ydata1 = listtempA2[5] #y2
        ydata2 = listtempA2[2] #x1
        ydata3 = listtempA2[3] #y1 
        ydata4 = listtempA2[0] #x0
        ydata5 = listtempA2[1] #y0
    data = [data0,data1,data2,data3,data4,data5,ydata0,ydata1,ydata2,ydata3,ydata4,ydata5]    
    return data
    #return data0,data1,data2,data3,data4,data5,ydata0,ydata1,ydata2,ydata3,ydata4,ydata5

###求平移###
def delhanshu(data0,data1,ydata0,ydata1):
    delx = ydata0-data0
    dely = ydata1-data1
        
    print('delx = ', delx)
    print('dely = ', dely)

    return delx,dely

def pingyiliang(data):
    data0 = data[0]
    data1 = data[1]
    data2 = data[2]
    data3 = data[3]
    data4 = data[4]
    data5 = data[5]
    
    ydata0 = data[6]
    ydata1 = data[7]
    ydata2 = data[8]
    ydata3 = data[9]
    ydata4 = data[10]
    ydata5 = data[11]
    
    delx1,dely1 =  delhanshu(data0,data1,ydata0,ydata1)   
    delx2,dely2 =  delhanshu(data2,data3,ydata2,ydata3) 
    delx3,dely3 =  delhanshu(data4,data5,ydata4,ydata5)
    delx = round((delx1+delx2+delx3)/3)
    dely = round((dely1+dely2+dely3)/3)
    return delx,dely
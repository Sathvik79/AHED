# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 08:13:36 2020

@author: Sathvik
"""

# Mean Residence and E-Curve

import matplotlib.pyplot as plt
#import numpy as np
from scipy.integrate import simps
from numpy import trapz


t = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14]
c = [0, 1, 5, 8, 10, 8, 6, 4, 3, 2.2, 1.5, 0.6, 0]
totalNum = 0
totalDeno = 0
eCurve = []


def CCurve():
    global totalNum
    global totalDeno
    dtI = 5
    for i in range(len(t)):   # tI = (t[i]*c[i]*dtI)/(c[i]*dtI)
        tNum = t[i]*c[i]*dtI
        totalNum += tNum
        tDeno = c[i]*dtI
        totalDeno += tDeno
    tI = totalNum/totalDeno
    print('ti = ', tI)
    print('ci = ', totalDeno)


def ECurve():
    global eCurve
    for i in range(len(c)):
       e = c[i]/totalDeno
       eCurve.append(e)
    print(eCurve)
    

def graphCCurve():
    plt.plot(t, c)
    plt.xlim(0, 15)
    plt.ylim(0, 10)
    plt.show()


def graphECurve():
    plt.plot(t, eCurve)
    plt.xlim(0, 35)
    plt.ylim(0, 0.06)
    plt.xlabel('t (s)')
    plt.ylabel('E-Curve (s-1)')
    plt.title('E-Curve')
    plt.show()


def area():
    '''
    totalArea = 0
    eCurve1 = eCurve
    eCurve2 = eCurve[1:]
    t1 = t
    t2 = t[1:len(t)]
    for i in range(len(eCurve)):
        total = (eCurve2[i]-eCurve1[i])*((t2[i]-t1[i])/2)
        totalArea += total
    print(totalArea)'''
    print('Area', trapz(c, dx = 35))
    print('Area', simps(c, dx = 35))
    print('Area', trapz(eCurve, dx = 35))
    print('Area', simps(eCurve, dx = 35))

    
CCurve()
ECurve()
graphCCurve()
graphECurve()
area()
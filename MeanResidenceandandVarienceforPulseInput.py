# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 21:40:45 2020

@author: Sathvik
"""

# Mean Residence Time and Varience by Pulse Input

import matplotlib.pyplot as plt
from scipy.integrate import simps
from numpy import trapz

t = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14]
E = [0, 0.02, 0.10, 0.16, 0.20, 0.16, 0.12, 0.08, 0.06, 0.044, 0.03, 0.012, 0]
tE = []
tEDashTotal = 0

def tEMulti():
    global tEDashTotal
    for i in range(len(t)):
        tETotal = t[i] * E[i]
        tEDashTotal += tETotal
        tE.append(tETotal)
    print(tE)
    print(tEDashTotal)
    
    
def area():
    plt.plot(t, E)
    plt.xlim(0, 15)
    plt.ylim(0, 0.20)
    plt.show()
    print('Area', trapz(t, dx = 15))
    print('Area', simps(t, dx = 14))
    print('Area', trapz(E, dx = 0.2))
    print('Area', simps(E, dx = 0.2))
    
    
tEMulti()
area()
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:45:53 2020

@author: rsanchez
"""

import numpy as np
import matplotlib.pyplot as plt
from context import camis
import time

h = 1.0

#x = np.array([0.5,np.sqrt(3)/2])*h
#xj = np.array([0.0, 0.0])*h
#xk = np.array([1.0, 0.0])*h

x = np.array([1.0,np.sqrt(3)/2])*h
xj = np.array([0.0, np.sqrt(3)/2])*h
xk = np.array([0.5, 0.0])*h

Tj = 1.0
Tk = 1.0
Q1 = 1.0
Q2 = 1.0
D1 = 0.0
D2 = 0.0
aspect = np.array([1.0, 0.0])


start = time.time()
for i in range(100):
    optimizedT,dirT = camis.optimizeCost(x,xj,xk,Tj,Tk,i,i,D1,D2,aspect)
print(time.time() - start)

start = time.time()
for i in range(100):
    computedT,computedDir = camis.getEikonalCost(x,xj,xk,Tj,Tk,i)
print(time.time() - start)
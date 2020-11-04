# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 18:53:46 2020

@author: Richi
"""

import numpy as np
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# PLANE FITTING FROM GRID POINTS

XX,YY = np.meshgrid(np.linspace(-0.5,0.5,11), np.linspace(-0.5,0.5,11))

ZZ = (XX+0.5)**2+(YY+0.5)**2

Xarray = XX.flatten()
Yarray = YY.flatten()
Zarray = ZZ.flatten()

A = np.c_[Xarray, Yarray, np.ones(Xarray.size)]
C,_,_,_ = scipy.linalg.lstsq(A, Zarray)

Zplane = C[0]*XX + C[1]*YY + C[2]

fig = plt.figure()
ax = fig.gca(projection='3d')
#ax.plot_surface(Xarray, Yarray, Zplane, rstride=1, cstride=1, alpha=0.2)
ax.scatter(Xarray, Yarray, Zarray, c='r', s=50)
ax.scatter(Xarray, Yarray, Zplane, c='b', s=50)
plt.xlabel('X')
plt.ylabel('Y')
ax.set_zlabel('Z')
ax.axis('equal')
ax.axis('tight')
plt.show()
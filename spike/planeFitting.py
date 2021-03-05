"""
MIT License
-----------

Copyright (c) 2021 University of Malaga
Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

Authors: J. Ricardo Sánchez Ibáñez, Carlos J. Pérez del Pulgar
Affiliation: Department of Systems Engineering and Automation
Space Robotics Lab (www.uma.es/space-robotics)
"""
# -*- coding: utf-8 -*-

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
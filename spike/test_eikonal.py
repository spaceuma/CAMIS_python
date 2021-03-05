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
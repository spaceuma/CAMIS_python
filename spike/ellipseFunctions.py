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
import scipy

deg2rad = np.pi/180
rad2deg = 180/np.pi

def getInverseProjectedSpeed(steepness_deg, beta_deg):
    steepness = deg2rad*steepness_deg
    beta = deg2rad*beta_deg
    return np.sqrt(np.cos(beta)**2 / np.cos(steepness)**2 + np.sin(beta)**2)

def getSlipFactor(srA,srB,saA,saB,steepness_deg,beta_deg):
    beta = deg2rad*beta_deg
    return np.sqrt(np.cos(beta)**2 / (1 - srA*np.exp(srB*steepness_deg))**2 + \
                   np.sin(beta)**2 / np.cos(saA*np.exp(saB*steepness_deg)*deg2rad)**2)

alpha = 25.0
beta_deg_array = np.arange(0, 360, 1.0)

pISpeed = []
slipFactor = []
cost = []


for beta_deg in beta_deg_array:
    pISpeed.append(getInverseProjectedSpeed(alpha,beta_deg))
    slipFactor.append(getSlipFactor(0.07,0.10,1.32,0.16,alpha,beta_deg))

for i in range(len(slipFactor)):
    cost.append(1 / (pISpeed[i]*slipFactor[i]))

fig = plt.figure()
axes5 = plt.subplot(111, projection='polar')
axes5.set_facecolor('xkcd:light blue')
axes5.plot(beta_deg_array*deg2rad, cost, 'xkcd:sea blue', lw = 2)
fig.tight_layout()
plt.show()
        
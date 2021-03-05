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

deg2rad = np.pi/180
rad2deg = 180/np.pi



steepnessArray = np.linspace(0,45,90+2)*deg2rad

rho = .4
steepnessThreshold = np.arctan(rho)

Rd = np.linspace(2,10,5)

descentTan = np.tan(steepnessArray)/rho

descentLine = steepnessArray/rho


def bezierBrakingCost(alpha, rhoValue):
    alphaBreaking1 = np.arctan(rhoValue)
    alphaBreaking2 = np.arctan(2*rhoValue)
    if alpha > alphaBreaking2:
        return np.abs(1 - np.tan(alpha)/rho)
    else:
        if alphaBreaking1 == alphaBreaking2:
            t = alpha/alphaBreaking2
        else:
            t = (-2*alphaBreaking1 + np.sqrt(4*alphaBreaking1**2 + 4*alpha*(alphaBreaking2 - 2*alphaBreaking1)))/(2*(alphaBreaking2 - 2*alphaBreaking1))
#            t2 = (-2*alphaBreaking1 - np.sqrt(4*alphaBreaking1**2 + 4*alpha*(alphaBreaking2 - 2*alphaBreaking1)))/(2*(alphaBreaking2 - 2*alphaBreaking1))
#            if t1 >= 0.0 and t1 <= 1.0:
#                t = t1
#            if t2 >= 0.0 and t2 <= 1.0:
#                t = t2
        Cost = 1 - 2*t + 2*t**2
        return Cost


fig, ax = plt.subplots(figsize=(4, 4),constrained_layout=True)
brakingFunction1 = np.zeros_like(steepnessArray)
ascentFunction = 1 + steepnessArray/rho
for i,steepness in enumerate(steepnessArray):
    brakingFunction1[i] = bezierBrakingCost(steepnessArray[i], rho)
ax.plot(steepnessArray*rad2deg, brakingFunction1)
ax.plot(steepnessArray*rad2deg, ascentFunction)
ax.plot(steepnessArray*rad2deg, np.abs(1 - descentTan), linestyle = 'dashed')
#ax.plot(steepnessArray*rad2deg, np.abs(1 - descentLine), linestyle = 'dashed')
ax.set_ylim([0.0, 5.0])
ax.set_xlabel('Steepness [degrees]')
ax.set_ylabel('Value')
ax.legend(('Rd = 2','Rd = 4','Rd = 6','Rd = 8','Rd = 10','tan/rho','steepness/rho'))
plt.show()



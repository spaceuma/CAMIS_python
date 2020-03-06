# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 15:27:50 2020

@author: rsanchez
"""

import numpy as np
import matplotlib.pyplot as plt

deg2rad = np.pi/180
rad2deg = 180/np.pi



steepnessArray = np.linspace(0,45,90.0+2)*deg2rad

steepnessMaxRisk = 30.0*deg2rad
steepnessMinRisk = 25.0*deg2rad
riskGain = 5.0

Kr = 0.3
Ko = Kr

rho = .7
steepnessBrake = np.arctan(rho)

#The intersection between risk function and ascent function is approximated
intersectionX = (steepnessMinRisk + (steepnessMaxRisk-steepnessMinRisk)/riskGain)/(1-(steepnessMaxRisk-steepnessMinRisk)/(rho*riskGain))
intersectionY = 1 + intersectionX/rho
intersectionPoint = np.array([intersectionX,intersectionY])

#The intersection between risk function and descent function is approximated
intersection2X = (steepnessMinRisk + (steepnessMaxRisk-steepnessMinRisk)/riskGain)/(1+(steepnessMaxRisk-steepnessMinRisk)/(rho*riskGain))
intersection2Y = 1 - intersectionX/rho
intersection2Point = np.array([intersection2X,intersection2Y])

initialPoint = np.array([(1-Ko)*intersectionX,1 + (1-Ko)*intersectionX/rho])
initial2Point = np.array([(1-Ko)*intersection2X,1 - (1-Ko)*intersection2X/rho])



riskPoint = np.array([intersectionX + Kr*(steepnessMaxRisk - intersectionX),intersectionY + Kr*(riskGain-intersectionY)])
riskPoint2 = np.array([intersection2X + Kr*(steepnessMaxRisk - intersection2X),intersection2Y + Kr*(riskGain-intersection2Y)])
#steepnessRisk = steepnessBrake + 5.0*deg2rad
#steepnessMin = min(steepnessRisk - 5.0*deg2rad, steepnessBrake)



Rd = np.linspace(2,10,5)

descentFunction = np.abs(1 - np.tan(steepnessArray)/rho)

descentLine = steepnessArray/rho


def bezierCost(x, P0, P1, P2):
    X1 = P0[0] - P1[0]
    X2 = P2[0] - P1[0]
    Y1 = P0[1] - P1[1]
    Y2 = P2[1] - P1[1]
    
    if x < P0[0]:
        return P0[1] + (x - P0[0])/(P1[0] - P0[0])*(P1[1] - P0[1])
    
    if x > P2[0]:
        return P2[1] + (x - P2[0])/(P2[0] - P1[0])*(P2[1] - P1[1])
#        return np.abs(1 - np.tan(alpha)/rho)
    else:
        if np.abs(X1 - X2) < 0.001:
            t = (x - P0[0] - P1[0] + X1)/(2*X1)
        else:
            t = (-2*X1 - np.sqrt(4*X1**2 - 4*(x-P1[0]-X1)*(-X1-X2)))/(2*(-X1-X2))
#            t2 = (-2*alphaBreaking1 - np.sqrt(4*alphaBreaking1**2 + 4*alpha*(alphaBreaking2 - 2*alphaBreaking1)))/(2*(alphaBreaking2 - 2*alphaBreaking1))
#            if t1 >= 0.0 and t1 <= 1.0:
#                t = t1
#            if t2 >= 0.0 and t2 <= 1.0:
#                t = t2
        Cost = P1[1] + (1 - t)**2*Y1 + t**2*Y2
        return Cost


costRisk = 2


fig, ax = plt.subplots(figsize=(4, 4),constrained_layout=True)
ascentCAMIS = np.zeros_like(steepnessArray)
descentCAMIS = np.zeros_like(steepnessArray)
ascentFunction = 1 + np.tan(steepnessArray)/rho
ascentApproxFunction = 1 + steepnessArray/rho

lateralFunction = np.sqrt(ascentFunction*descentFunction)

pitchRiskFunction = (steepnessArray - steepnessMinRisk)/(steepnessMaxRisk - steepnessMinRisk)*riskGain
rollRiskFunction = (steepnessArray - 20.0*deg2rad)/(5.0*deg2rad)*riskGain

ax.plot(steepnessArray*rad2deg, ascentFunction, linestyle = 'dashed', color = 'b')
ax.plot(steepnessArray*rad2deg, ascentApproxFunction, linestyle = 'dotted', color = 'b')
for i,steepness in enumerate(steepnessArray):
    ascentCAMIS[i] = bezierCost(steepnessArray[i], initialPoint, intersectionPoint, riskPoint)
    descentCAMIS[i] = bezierCost(steepnessArray[i], initial2Point, intersection2Point, riskPoint2)
ax.plot(steepnessArray*rad2deg, ascentCAMIS, color = 'b')
ax.plot(steepnessArray*rad2deg, descentFunction, linestyle = 'dashed', color = 'g')
ax.plot(steepnessArray*rad2deg, descentCAMIS, color = 'g')
ax.plot(steepnessArray*rad2deg, lateralFunction, linestyle = 'dashed')
ax.plot(steepnessArray*rad2deg, pitchRiskFunction, linestyle = 'solid')
ax.plot(steepnessArray*rad2deg, rollRiskFunction, linestyle = 'solid')
ax.plot(intersectionX*rad2deg, intersectionY, '*')
ax.plot(riskPoint[0]*rad2deg, riskPoint[1], '*')
ax.plot(initialPoint[0]*rad2deg, initialPoint[1], '*')
ax.plot(initial2Point[0]*rad2deg, initial2Point[1], '*')
ax.plot(intersection2Point[0]*rad2deg, intersection2Point[1], '*')
#ax.plot(steepnessArray*rad2deg, np.abs(1 - descentLine), linestyle = 'dashed')
ax.set_ylim([0.0, riskGain])
ax.set_xlabel('Steepness [degrees]')
ax.set_ylabel('Value')
ax.legend(('Ascent','Ascent Approximation','CAMIS ascent','Descent','Lateral','Pitch Tip-Over Risk','Roll Tip-Over Risk','Intersection','Risk Point', 'initialPoint'))
plt.show()
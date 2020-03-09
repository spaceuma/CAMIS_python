# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 15:27:50 2020

@author: rsanchez
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import optimize
from scipy.optimize import least_squares

deg2rad = np.pi/180
rad2deg = 180/np.pi


#Configurable Parameters
rho = .4
pitchRiskMax = 30.0*deg2rad
pitchRiskMin = pitchRiskMax - 5.0*deg2rad
blockRiskMax = 25.0*deg2rad
blockRiskMin = blockRiskMax - 5.0*deg2rad
rollRiskMax = 20.0*deg2rad
rollRiskMin = rollRiskMax - 5.0*deg2rad
riskGain = 5.0
Kr = 0.4
Kr = np.maximum(Kr, 0.15)

#Derived Parameters
brakeSteepness = np.arctan(rho)
brakePoint = np.array([brakeSteepness,0])

steepnessArray = np.linspace(0.0,45.0,90+2)*deg2rad



#The intersection between risk function and ascent function is approximated

intersectionX = scipy.optimize.fsolve(lambda x : 
    (x - blockRiskMin)/(blockRiskMax-blockRiskMin)*riskGain - 
    np.tan(x) - rho,0)
intersectionY = rho + np.tan(intersectionX)
ascentIntersection = np.array([intersectionX,intersectionY])

#The intersection between risk function and descent function is approximated
if brakeSteepness >= pitchRiskMin:
    intersection2X = scipy.optimize.fsolve(lambda x : 
        (x - pitchRiskMin)/(pitchRiskMax-pitchRiskMin)*riskGain + 
        np.tan(x) - rho,0)
    intersection2Y = rho - np.tan(intersection2X)
else:
    intersection2X = scipy.optimize.fsolve(lambda x : 
        (x - pitchRiskMin)/(pitchRiskMax-pitchRiskMin)*riskGain + 
        -np.tan(x) + rho,0)
    intersection2Y = - rho + np.tan(intersection2X)
descentIntersection = np.array([intersection2X,intersection2Y])

#The intersection between risk function and lateral function is approximated
if brakeSteepness >= rollRiskMin:
    intersection3X = scipy.optimize.fsolve(lambda x : 
        (x - rollRiskMin)/(rollRiskMax-rollRiskMin)*riskGain - 
        np.sqrt(rho**2  - np.tan(x)**2),0)
    intersection3Y = np.sqrt(rho**2  - np.tan(intersection3X)**2)
else:
    intersection3X = scipy.optimize.fsolve(lambda x : 
        (x - rollRiskMin)/(rollRiskMax-rollRiskMin)*riskGain - 
        np.sqrt(np.tan(x)**2 - rho**2 ),rollRiskMax)
    intersection3Y = np.sqrt(np.tan(intersection3X)**2 - rho**2 )
lateralIntersection = np.array([intersection3X,intersection3Y])

#The first Bezier points
initialPoint = np.array([(1-Kr)*intersectionX,
                         rho + np.tan((1-Kr)*intersectionX)])

if brakeSteepness >= pitchRiskMin:
    initial2Point = np.array([(1-Kr)*intersection2X,
                              rho - np.tan((1-Kr)*intersection2X)])
else:
    initial2Point = np.array([(1-Kr)*brakeSteepness,
                              rho - np.tan((1-Kr)*brakeSteepness)])
    
if brakeSteepness >= rollRiskMin:
    initial3Point = np.array([(1-Kr)*intersection3X,
                              np.sqrt(rho**2  - np.tan((1-Kr)*intersection3X)**2)])
else:
    initial3Point = np.array([(1-Kr)*brakeSteepness, 
                              np.sqrt(rho**2  - np.tan((1-Kr)*brakeSteepness)**2)])

# The last Bezier points
riskPoint = np.array([intersectionX + Kr*(blockRiskMax - intersectionX),
                      intersectionY + Kr*(riskGain-intersectionY)])
riskPoint2 = np.array([intersection2X + Kr*(pitchRiskMax - intersection2X),
                       intersection2Y + Kr*(riskGain-intersection2Y)])
riskPoint3 = np.array([intersection3X + Kr*(rollRiskMax - intersection3X),
                       intersection3Y + Kr*(riskGain-intersection3Y)])

descentFunction = np.abs(rho - np.tan(steepnessArray))
ascentFunction = rho + np.tan(steepnessArray)


def bezierCost(x, P0, P1, P2):
    X1 = P0[0] - P1[0]
    X2 = P2[0] - P1[0]
    Y1 = P0[1] - P1[1]
    Y2 = P2[1] - P1[1]
    
    if x < P0[0]:
        return P0[1] - (P0[0] - x)/(P1[0] - P0[0])*(P1[1] - P0[1])
    
    if x > P2[0]:
        return P2[1] + (x - P2[0])/(P2[0] - P1[0])*(P2[1] - P1[1])
    else:
        if np.abs(X1 - X2) < 0.001:
            t = (x - P0[0] - P1[0] + X1)/(2*X1)
        else:
            t = (-2*X1 - np.sqrt(4*X1**2 - 4*(x-P1[0]-X1)*(-X1-X2)))/(2*(-X1-X2))
        Cost = P1[1] + (1 - t)**2*Y1 + t**2*Y2
        return Cost
    
def bezierCubicCost(steepness,P0,P1,P2,P3):
#    t = scipy.optimize.fsolve(lambda x : 
#    (1-x)**3*P0[0] + 3*(1-x)**2*x*P1[0] + 3*(1-x)*x**2*P2[0]+x**3*P3[0],0)
    res = scipy.optimize.least_squares(lambda x : 
    (1-x)**3*P0[0] + 3*(1-x)**2*x*P1[0] + 3*(1-x)*x**2*P2[0]+x**3*P3[0] - steepness,0.5, bounds = (0,1))
    t = res.x
    return (1-t)**3*P0[1] + 3*(1-t)**2*t*P1[1] + 3*(1-t)*t**2*P2[1]+t**3*P3[1]

plt.style.use('seaborn')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, ax = plt.subplots(figsize=(5, 4),constrained_layout=True, num = 'seaborn')
ascentCAMIS = np.zeros_like(steepnessArray)
descentCAMIS = np.zeros_like(steepnessArray)
lateralCAMIS = np.zeros_like(steepnessArray)


lateralFunction = np.sqrt(ascentFunction*descentFunction)

pitchRiskFunction = (steepnessArray - pitchRiskMin)/(pitchRiskMax - pitchRiskMin)*riskGain
rollRiskFunction = (steepnessArray - rollRiskMin)/(rollRiskMax - rollRiskMin)*riskGain
blockRiskFunction = (steepnessArray - blockRiskMin)/(blockRiskMax - blockRiskMin)*riskGain



ax.plot(steepnessArray*rad2deg, ascentFunction, linestyle = 'dashed', color = 'b')
for i,steepness in enumerate(steepnessArray):
    
    if steepness < initialPoint[0]:
        ascentCAMIS[i] = rho + np.tan(steepnessArray[i])
    else:
        if steepness > riskPoint[0]:
            ascentCAMIS[i] = (steepness - blockRiskMin)/(blockRiskMax - blockRiskMin)*riskGain
        else:
            ascentCAMIS[i] = bezierCost(steepnessArray[i], initialPoint, ascentIntersection, riskPoint)
        
    if steepness < initial2Point[0]:
        descentCAMIS[i] = np.abs(rho - np.tan(steepnessArray[i])) #TODO: take care of braking!!
    else:
        if steepness >= riskPoint2[0]:
            descentCAMIS[i] = (steepness - pitchRiskMin)/(pitchRiskMax - pitchRiskMin)*riskGain
        else:
            if brakeSteepness >= pitchRiskMin:
                descentCAMIS[i] = bezierCost(steepnessArray[i], initial2Point, descentIntersection, riskPoint2)
            else:
                descentCAMIS[i] = bezierCubicCost(steepnessArray[i], initial2Point, brakePoint,descentIntersection, riskPoint2)
                
    if steepness < initial3Point[0]:
        lateralCAMIS[i] = np.sqrt(rho**2  - np.tan(steepnessArray[i])**2)
    else:
        if steepness >= riskPoint3[0]:
            lateralCAMIS[i] = (steepness - rollRiskMin)/(rollRiskMax - rollRiskMin)*riskGain
        else:
            if brakeSteepness >= rollRiskMin:
                lateralCAMIS[i] = bezierCost(steepnessArray[i], initial3Point, lateralIntersection, riskPoint3)
            else:
                lateralCAMIS[i] = bezierCost(steepnessArray[i], initial3Point, lateralIntersection, riskPoint3)
                
ax.plot(steepnessArray*rad2deg, ascentCAMIS, color = 'b')
ax.plot(intersectionX*rad2deg, intersectionY, '*', color = 'b')

ax.plot(steepnessArray*rad2deg, pitchRiskFunction, linestyle = 'dotted', color = 'r')
ax.plot(steepnessArray*rad2deg, rollRiskFunction, linestyle = 'dotted', color = 'm')
ax.plot(steepnessArray*rad2deg, blockRiskFunction, linestyle = 'dotted', color = 'y')
ax.plot(riskPoint[0]*rad2deg, riskPoint[1], '*', color = 'b')


ax.plot(initialPoint[0]*rad2deg, initialPoint[1], '*', color = 'b')
ax.plot(lateralIntersection[0]*rad2deg, lateralIntersection[1], '*', color = 'orange')
ax.set_ylim([0.0, riskGain])
ax.set_xlabel('Steepness [degrees]')
ax.set_ylabel('Value')
ax.legend(('Ascent','Ascent Approximation','CAMIS ascent','Descent','Lateral','Pitch Tip-Over Risk','Roll Tip-Over Risk','Intersection','Risk Point', 'initialPoint'))
plt.show()

anisoAD = ascentCAMIS/descentCAMIS
anisoLD = lateralCAMIS/descentCAMIS

fig1, ax1 = plt.subplots(figsize=(5, 4),constrained_layout=True)
ax1.plot(steepnessArray*rad2deg, ascentFunction, linestyle = 'dashed', color = 'b')
ax1.plot(steepnessArray*rad2deg, descentFunction, linestyle = 'dashed', color = 'g')
ax1.plot(steepnessArray*rad2deg, lateralFunction, linestyle = 'dashed', color = 'orange')
ax1.plot(steepnessArray*rad2deg, ascentCAMIS, color = 'b')
ax1.plot(steepnessArray*rad2deg, descentCAMIS, color = 'g')
ax1.plot(steepnessArray*rad2deg, lateralCAMIS, color = 'orange')
ax1.plot(steepnessArray*rad2deg, blockRiskFunction, linestyle = 'dotted', color = 'y')
ax1.plot(steepnessArray*rad2deg, pitchRiskFunction, linestyle = 'dotted', color = 'r')
ax1.plot(steepnessArray*rad2deg, rollRiskFunction, linestyle = 'dotted', color = 'm')
ax1.plot(initialPoint[0]*rad2deg, initialPoint[1], 'o', color = 'b')
ax1.plot(initial2Point[0]*rad2deg, initial2Point[1], 'o', color = 'g')
ax1.plot(initial3Point[0]*rad2deg, initial3Point[1], 'o', color = 'orange')
ax1.plot(intersectionX*rad2deg, intersectionY, 'o', color = 'b')
ax1.plot(descentIntersection[0]*rad2deg, descentIntersection[1], 'o', color = 'g')
ax1.plot(lateralIntersection[0]*rad2deg, lateralIntersection[1], 'o', color = 'orange')
ax1.plot(riskPoint3[0]*rad2deg, riskPoint3[1], 'o', color = 'orange')
ax1.plot(riskPoint[0]*rad2deg, riskPoint[1], 'o', color = 'b')
ax1.plot(brakePoint[0]*rad2deg, brakePoint[1], 'o', color = 'g')
ax1.plot(riskPoint2[0]*rad2deg, riskPoint2[1], 'o', color = 'g')
ax1.plot(riskPoint3[0]*rad2deg, riskPoint3[1], 'o', color = 'orange')
ax1.legend(('ρ + tan α','|ρ - tan α|','((ρ + tan α)|ρ - tan α|)^.5','$R_a$','$R_d$','$R_l$','Block Risk Line', 
            'Pitch Risk Line','Roll Risk Line','Ascent Bezier points', 'Descent Bezier points', 'Lateral Bezier points'))
ax1.set_xlim([0.0, 30.0])
ax1.set_ylim([0.0, 3.0])
ax1.set_xlabel('Steepness [degrees]')
ax1.set_ylabel('R')

fig2, ax2 = plt.subplots(figsize=(4, 4),constrained_layout=True)
ax2.plot(steepnessArray*rad2deg, anisoAD, linestyle = 'solid')
ax2.plot(steepnessArray*rad2deg, anisoLD, linestyle = 'solid')
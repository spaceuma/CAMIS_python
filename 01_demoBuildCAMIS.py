# -*- coding: utf-8 -*-
"""
DEMO script: building CAMIS upon experimental data
@author: J.Ricardo Sanchez Ibanez (ricardosan@uma.es), Github: @Ryk-San

This script is a reference demo to understand how CAMIS is built using data
obtained by a mobile robot.

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import lib.cuadriga_reader as cr
import lib.camis as camis


"01 - The CSV files collected during preliminar experiments are loaded"
csvFiles = ['cuadrigaData/20190531/JornadasRescate01.txt',
            'cuadrigaData/20190531/JornadasRescate02.txt',
            'cuadrigaData/20190531/JornadasRescate03.txt',
            'cuadrigaData/201811/slopeTest01.txt',
            'cuadrigaData/201811/slopeTest02.txt',
            'cuadrigaData/201811/slopeTest04.txt',
            'cuadrigaData/201811/slopeTest05.txt',
            'cuadrigaData/201811/slopeTest06.txt',
            'cuadrigaData/201811/slopeTest07.txt',
            'cuadrigaData/201811/slopeTest08.txt',
            'cuadrigaData/201811/slopeTest09.txt',
            'cuadrigaData/201811/slopeTest10.txt'
            ]


"02 - CAMIS input data is extracted"
betaX = []
betaY = []
cost = []
gradient = []
for file in csvFiles:
    " Cost C according to robot orientation Roll-Pitch-Yaw "
    Roll, Pitch, Yaw, C = cr.readCuadrigaData(file)
    " Conversion from Roll-Pitch-Yaw to Gradient-Beta_angle "
    G, Bx, By = camis.rpy2ab(Roll,Pitch,Yaw)
    betaX = betaX + Bx
    betaY = betaY + By
    cost = cost + C
    gradient = gradient + G
    
    
"03 - The functions Cd, Ca, Cl1 and Cl2 are obtained through a fitting process"
descentCurrent = []
descentGradient = []
ascentCurrent = []
ascentGradient = []
lateral1Current = []
lateral1Gradient = []
lateral2Current = []
lateral2Gradient = []
plainCurrent = []
for i,b in enumerate(betaX):
    if betaX[i] > 0:
        descentCurrent.append(cost[i])
        descentGradient.append(gradient[i])
    else:
        ascentCurrent.append(cost[i])
        ascentGradient.append(gradient[i])
    if betaY[i] > 0:
        lateral1Current.append(cost[i])
        lateral1Gradient.append(gradient[i])
    else:
        lateral2Current.append(cost[i])
        lateral2Gradient.append(gradient[i])
    if gradient[i] < 7:
        plainCurrent.append(cost[i])
 
plainCost = np.mean(plainCurrent)

def costFunction(magnitude, a, b):
    return a * magnitude**2 + b*magnitude + plainCost

dpopt, dpcov = curve_fit(costFunction, descentGradient, descentCurrent)
apopt, apcov = curve_fit(costFunction, ascentGradient, ascentCurrent)
l1popt, l1pcov = curve_fit(costFunction, lateral1Gradient, lateral1Current)
l2popt, l2pcov = curve_fit(costFunction, lateral2Gradient, lateral2Current)

gradient = np.asarray(gradient)
linearGradient = np.linspace(0,30,31)

cm = plt.get_cmap("RdYlGn")
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
xs = np.multiply(betaX,cost)
ys = np.multiply(betaY,cost)
col = np.zeros_like(xs)
for i,s in enumerate(xs):
    col[i] = np.sqrt(np.power(xs[i],2)+np.power(ys[i],2))
ax.scatter(xs, ys, gradient,c = col, alpha = .1)
ax.plot(costFunction(linearGradient, dpopt[0], dpopt[1]), np.zeros_like(linearGradient), linearGradient, color = 'r')
ax.plot(-costFunction(linearGradient, apopt[0], apopt[1]), np.zeros_like(linearGradient), linearGradient, color = 'r')
ax.plot(np.zeros_like(linearGradient), costFunction(linearGradient, l1popt[0], l1popt[1]), linearGradient, color = 'r')
ax.plot(np.zeros_like(linearGradient), -costFunction(linearGradient, l2popt[0], l2popt[1]), linearGradient, color = 'r')
ax.set_aspect('equal')

heading = np.arange(0, 2*np.pi, 0.01)
aspect = 0  

Xs = []
Ys = []

for g in linearGradient:
    for theta in heading:
        Cd = costFunction(g, dpopt[0], dpopt[1])
        Ca = costFunction(g, apopt[0], apopt[1])
        Cl1 = costFunction(g, l1popt[0], l1popt[1])
        Cl2 = costFunction(g, l2popt[0], l2popt[1])
        preCost,T = camis.computeCAMIScost(theta,aspect,Cd,Ca,Cl1,Cl2)
        Xs.append(T[0]*preCost)
        Ys.append(T[1]*preCost)
    ax.plot(Xs, Ys, g, color=plt.cm.jet(float(g)/25))
    ax.plot(Xs, Ys, 0, color = 'k')
    Xs = []
    Ys = []


plt.show()



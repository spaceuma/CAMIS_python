# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 11:43:38 2020

@author: Richi
"""
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
        
# -*- coding: utf-8 -*-
"""
DEMO script: building CAMIS based on Cuadriga experimental data
@author: J.Ricardo Sanchez Ibanez (ricardosan@uma.es), Github: @Ryk-San

This script is a reference demo to understand how CAMIS is built using data
obtained by a mobile robot.

"""
import numpy as np
import matplotlib.pyplot as plt
import lib.camis as cm
import lib.anisotropic_path_planning as ap

XMap = np.loadtxt(open('UMARescueArea_VCM_X.csv',\
                                     "rb"), delimiter=" ", skiprows=0)
YMap = np.loadtxt(open('UMARescueArea_VCM_Y.csv',\
                                     "rb"), delimiter=" ", skiprows=0)
AnisotropyMap = np.loadtxt(open('UMARescueArea_VCM_Anisotropy.csv',\
                                     "rb"), delimiter=" ", skiprows=0)

VCMap = np.zeros([4,AnisotropyMap.shape[0],AnisotropyMap.shape[1]])

VCMap[0] = np.loadtxt(open('UMARescueArea_VCM_Q1.csv',\
                                     "rb"), delimiter=" ", skiprows=0)
VCMap[1] = np.loadtxt(open('UMARescueArea_VCM_Q2.csv',\
                                     "rb"), delimiter=" ", skiprows=0)
VCMap[2] = np.loadtxt(open('UMARescueArea_VCM_D1.csv',\
                                     "rb"), delimiter=" ", skiprows=0)
VCMap[3] = np.loadtxt(open('UMARescueArea_VCM_D2.csv',\
                                     "rb"), delimiter=" ", skiprows=0)

AspectMap = np.zeros([2,AnisotropyMap.shape[0],AnisotropyMap.shape[1]])

AspectMap[0] = np.loadtxt(open('UMARescueArea_VCM_AspectX.csv',\
                                     "rb"), delimiter=" ", skiprows=0)

AspectMap[1] = np.loadtxt(open('UMARescueArea_VCM_AspectY.csv',\
                                     "rb"), delimiter=" ", skiprows=0)


#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(VCMap[0], 100)
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()
#
#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(XMap, YMap, np.abs(np.arctan2(AspectMap[1],AspectMap[0])), 100, cmap = 'hsv')
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()

start = np.asarray([190,75])
goal = np.asarray([400,25])

x = np.asarray([0,0])
xj = np.asarray([0,1])
xk = np.asarray([1,0])
#T,dire = ap.optimizeCost(x,xj,xk,5,5,4,4,0,0,0)

globalRes = XMap[0,1] - XMap[0,0]

TmapG,TmapS, dirMapG, dirMapS, nodeLink, stateMapG,stateMapS = ap.computeTmap(VCMap,AspectMap,AnisotropyMap,goal,start,XMap,YMap,globalRes)
#TmapG,TmapS = ap.computeTmap(VCMap,AspectMap,AnisotropyMap,goal,start)

fig, axes = plt.subplots(constrained_layout=True)
cc = axes.contour(XMap,YMap, TmapG, 50)
fig.colorbar(cc,location='bottom')
axes.set_aspect('equal')
plt.show()







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

from time import time

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

#AnisotropyMap[:] = 2
##
#VCMap[0] = 2*AnisotropyMap
#VCMap[1] = AnisotropyMap
#VCMap[2] = 0*AnisotropyMap
#VCMap[3] = 0*AnisotropyMap
#AspectMap[0] = AnisotropyMap
#AspectMap[1] = 0*AnisotropyMap
#
#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(XMap, YMap, AnisotropyMap, 100, cmap = 'plasma')
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()



#VCMap[1] = VCMap[0]
#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(np.sqrt(VCMap[1])-np.sqrt(VCMap[0]), 100, cmap = 'plasma')
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()
#
#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(XMap, YMap, AnisotropyMap, 100, cmap = 'plasma')
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()

#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(XMap, YMap, np.abs(np.arctan2(AspectMap[1],AspectMap[0])), 100, cmap = 'hsv')
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()

goal = np.asarray([200,300])
start = np.asarray([140,165])
#
#x = np.asarray([0,0])
#xj = np.asarray([0,1])
#xk = np.asarray([1,0])
##T,dire = ap.optimizeCost(x,xj,xk,5,5,4,4,0,0,0)

globalRes = XMap[0,1] - XMap[0,0]
init = time()
TmapG,TmapS, dirMapG, dirMapS, nodeLink, stateMapG,stateMapS = ap.computeTmap(VCMap,AspectMap,AnisotropyMap,goal,start,XMap,YMap,globalRes)
print('Elapsed time to compute the Total Cost Map: '+str(time()-init))
###TmapG,TmapS = ap.computeTmap(VCMap,AspectMap,AnisotropyMap,goal,start)
##
fig, axes = plt.subplots(constrained_layout=True)
cc = axes.contourf(XMap,YMap, stateMapG, 50)
fig.colorbar(cc,location='bottom')
axes.set_aspect('equal')
plt.show()

Tmap = np.zeros_like(TmapG)
Tmap[:] = TmapS
for j in range(Tmap.shape[0]):
    for i in range(Tmap.shape[1]):
        if np.isinf(TmapG[j,i]):
            pass
        else:
            Tmap[j,i] = TmapS[nodeLink[1],nodeLink[0]] + TmapG[nodeLink[1],nodeLink[0]] - TmapG[j,i]

fig, ax = plt.subplots()
ax.contourf(XMap, YMap, Tmap, 100, cmap = 'magma', alpha = .5)
ax.contour(XMap, YMap, Tmap, 100, cmap = 'magma')
ax.set_aspect('equal')
plt.show()


#
fig, axes = plt.subplots(constrained_layout=True)
cc = axes.contourf(XMap, YMap, dirMapG, 100, cmap = 'hsv')
fig.colorbar(cc,location='bottom')
axes.set_aspect('equal')
plt.show()





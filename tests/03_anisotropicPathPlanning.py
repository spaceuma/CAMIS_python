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
XY2I = np.loadtxt(open('UMARescueArea_VCM_XY2I.csv',\
                                     "rb"), delimiter=" ", skiprows=0)
XY2J = np.loadtxt(open('UMARescueArea_VCM_XY2J.csv',\
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

utmOrigin = np.loadtxt(open('UMARescueArea_VCM_utmOrigin.csv',\
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
#cc = axes.contourf(AnisotropyMap, 100, cmap = 'plasma')
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

#CASE A
#start = np.asarray([60,60])
#goal = np.asarray([145,125])

#CASE B
#start = np.asarray([175,60])
#goal = np.asarray([45,215])

#CASE C
#goal = np.asarray([35,140])
#start = np.asarray([65,150])

start = np.asarray([100,50])
#start = np.asarray([30,175])
goal = np.asarray([100,100])

#
#x = np.asarray([0,0])
#xj = np.asarray([0,1])
#xk = np.asarray([1,0])
##T,dire = ap.optimizeCost(x,xj,xk,5,5,4,4,0,0,0)

globalRes = XMap[0,1] - XMap[0,0]
init = time()
Tmap1, dirMap1, stateMap1, maxAnisoMap1 = ap.computeTmap(VCMap,AspectMap,\
                                                     AnisotropyMap,goal,start,\
                                                     XMap,YMap,globalRes)

Tmap2, dirMap2, stateMap2, maxAnisoMap2 = ap.computeTmap(VCMap,AspectMap,\
                                                     AnisotropyMap,start,goal,\
                                                     XMap,YMap,globalRes)

print('Elapsed time to compute the Total Cost Map: '+str(time()-init))
###TmapG,TmapS = ap.computeTmap(VCMap,AspectMap,AnisotropyMap,goal,start)
##

IJ2XY = np.zeros([2,XMap.shape[0],YMap.shape[1]])
IJ2XY[0] = XMap
IJ2XY[1] = YMap
XY2IJ = np.zeros([2,XY2I.shape[0],XY2J.shape[1]])
XY2IJ[0] = XY2I
XY2IJ[1] = XY2J


startWaypoint = IJ2XY[:,start[1],start[0]]
goalWaypoint = IJ2XY[:,goal[1],goal[0]]
path1,u1 = ap.getPath(dirMap1, IJ2XY, XY2IJ, startWaypoint, goalWaypoint, utmOrigin[0], utmOrigin[1], globalRes)
path1 = np.asarray(path1)
u1 = np.asarray(u1)

u1x = IJ2XY[0,u1[:,1],u1[:,0]]
u1y = IJ2XY[1,u1[:,1],u1[:,0]]

path2,u2 = ap.getPath(dirMap2, IJ2XY, XY2IJ, goalWaypoint, startWaypoint, utmOrigin[0], utmOrigin[1], globalRes)
path2 = np.asarray(path2)
u2 = np.asarray(u2)

u2x = IJ2XY[0,u2[:,1],u2[:,0]]
u2y = IJ2XY[1,u2[:,1],u2[:,0]]

fig, axes = plt.subplots(constrained_layout=True)
cc = axes.contourf(AnisotropyMap, 50)
#axes.plot(path1[:,0],path1[:,1],'r')
#axes.plot(path2[:,0],path2[:,1],'g')
fig.colorbar(cc,location='bottom')
axes.set_aspect('equal')
plt.show()

fig, axes = plt.subplots(constrained_layout=True)
cc = axes.contourf(XMap, YMap, np.sqrt(VCMap[0]), 50)
axes.plot(path1[:,0],path1[:,1],'r')
axes.plot(path2[:,0],path2[:,1],'g')
fig.colorbar(cc,location='bottom')
axes.set_aspect('equal')
plt.show()

fig, axes = plt.subplots(constrained_layout=True)
cc = axes.contourf(XMap, YMap, VCMap[2], 50)
axes.plot(path1[:,0],path1[:,1],'r')
axes.plot(path2[:,0],path2[:,1],'g')
fig.colorbar(cc,location='bottom')
axes.set_aspect('equal')
plt.show()


fig, ax = plt.subplots(constrained_layout=True)
cc = ax.contourf(XMap, YMap, Tmap1, np.linspace(0,1000), cmap = 'magma', alpha = .5, vmin=0, vmax=600)
#ax.contour(XMap, YMap, Tmap1, 100, cmap = 'magma')
ax.plot(path1[:,0],path1[:,1],'r')
fig.colorbar(cc,location='bottom')
ax.set_aspect('equal')
plt.show()

fig, ax = plt.subplots(constrained_layout=True)
cc = ax.contourf(XMap, YMap, Tmap2, 100, cmap = 'magma', alpha = .5)
ax.contour(XMap, YMap, Tmap2, 100, cmap = 'magma')
#ax.plot(path2[:,0],path2[:,1],'g')
fig.colorbar(cc,location='bottom')
ax.set_aspect('equal')
plt.show()


#
#
##
#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(XMap, YMap, dirMap, 100, cmap = 'hsv')
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()
#
#
fig, axes = plt.subplots()
cc = axes.quiver(XMap, YMap, np.cos(dirMap1), np.sin(dirMap1),scale=20)
axes.scatter(u1x,u1y)
axes.scatter(path1[:,0],path1[:,1])
axes.plot(path1[:,0],path1[:,1],'r')
#fig.colorbar(cc,location='bottom')
axes.set_aspect('equal')
plt.show()


fig, axes = plt.subplots()
cc = axes.quiver(XMap, YMap, np.cos(dirMap2), np.sin(dirMap2),scale=20)
axes.scatter(u2x,u2y)
axes.scatter(path2[:,0],path2[:,1])
#axes.plot(path2[:,0],path2[:,1],'g')
#fig.colorbar(cc,location='bottom')
axes.set_aspect('equal')
plt.show()

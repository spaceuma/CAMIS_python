# -*- coding: utf-8 -*-
"""
DEMO script: building Vectorial Cost Map
@author: J.Ricardo Sanchez Ibanez (ricardosan@uma.es), Github: @Ryk-San

"""

import numpy as np
from time import time
from matplotlib import cbook
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt
from matplotlib.cbook import get_sample_data
import lib.camis as camis
import lib.hexgrid as hg
import lib.terrain_processing as tp
import math

from scipy import ndimage
import scipy.interpolate as interp

occupancyRadius = .5 #meters
globalRes = .5
slopeThreshold = 30

CdRoots,CaRoots,Cl1Roots,Cl2Roots,AniCoLUT = \
camis.readCamis('cuadriga_camis.csv')

hiRes_elevationMap = np.loadtxt(open("terrainData/UMARescueArea_10cmDEM.csv",\
                                     "rb"), delimiter=" ", skiprows=0)
hiRes_posX = np.loadtxt(open("terrainData/UMARescueArea_10cmPosX.csv",\
                                     "rb"), delimiter=" ", skiprows=0)
hiRes_posY = np.loadtxt(open("terrainData/UMARescueArea_10cmPosY.csv",\
                                     "rb"), delimiter=" ", skiprows=0)

hiRes = hiRes_posX[0,1] - hiRes_posX[0,0]

slopeMap, aspectMapX, aspectMapY, smoothDEM = tp.getConvSlopeMaps(hiRes_posX, hiRes_posY, hiRes_elevationMap,\
                                                       occupancyRadius, hiRes)
obstacleMap = slopeMap > slopeThreshold
obstacleMap = obstacleMap.astype(int)
obstacleMap[0,:] = 1
obstacleMap[-1,:] = 1
obstacleMap[:,0] = 1
obstacleMap[:,-1] = 1
proximityMap = hiRes*ndimage.morphology.distance_transform_edt(1-obstacleMap)
proximityMap = proximityMap - occupancyRadius
proximityMap[np.where(proximityMap[:]<0)] = 0


# Hexagonal Grid
triX, triY = hg.getHexGrid(hiRes_posX,hiRes_posY,globalRes)

points = np.zeros((hiRes_posX.size,2))

points[:,0] = hiRes_posX.flatten()
points[:,1] = hiRes_posY.flatten()

init = time()

tri_slopeMap = interp.griddata(points, slopeMap.flatten(), (triX, triY), method='nearest')
tri_slopeMap[np.where(triX < hiRes_posX[0,0])] = np.nan
tri_slopeMap[np.where(triX > hiRes_posX[-1,-1])] = np.nan
tri_slopeMap[np.where(triY < hiRes_posY[0,0])] = np.nan
tri_slopeMap[np.where(triY > hiRes_posY[-1,-1])] = np.nan

tri_VCMap = camis.getVectorialCostMap(\
                                      tri_slopeMap,CdRoots,CaRoots,Cl1Roots,\
                                      Cl2Roots,AniCoLUT)
#tri_elevationMap = interp.griddata(points, hiRes_elevationMap.flatten(), (triX-hiRes_posX[0,0], triY-hiRes_posY[0,0]), method='cubic')
tri_aspectMapX = interp.griddata(points, aspectMapX.flatten(), (triX, triY), method='nearest')
tri_aspectMapY = interp.griddata(points, aspectMapY.flatten(), (triX, triY), method='nearest')
tri_proximityMap = interp.griddata(points, proximityMap.flatten(), (triX, triY), method='nearest')

tri_aspectMapX[np.where(np.isnan(tri_slopeMap))] = np.nan
tri_aspectMapY[np.where(np.isnan(tri_slopeMap))] = np.nan
tri_proximityMap[np.where(np.isnan(tri_slopeMap))] = np.nan

obstacleMask = tri_proximityMap == 0

anisotropyMap = tri_VCMap[0][:][:]
anisotropyMap[obstacleMask] = np.inf
Q1 = tri_VCMap[1][:][:]
Q1[obstacleMask] = np.inf
Q2 = tri_VCMap[2][:][:]
Q2[obstacleMask] = np.inf
D1 = tri_VCMap[3][:][:]
D1[obstacleMask] = np.inf
D2 = tri_VCMap[4][:][:]
D2[obstacleMask] = np.inf

print('Elapsed time to compute all: '+str(time()-init))

np.savetxt('UMARescueArea_VCM_X.csv', triX, delimiter=" ")
np.savetxt('UMARescueArea_VCM_Y.csv', triY, delimiter=" ")
np.savetxt('UMARescueArea_VCM_Anisotropy.csv', anisotropyMap, delimiter=" ")
np.savetxt('UMARescueArea_VCM_Q1.csv', Q1, delimiter=" ")
np.savetxt('UMARescueArea_VCM_Q2.csv', Q2, delimiter=" ")
np.savetxt('UMARescueArea_VCM_D1.csv', D1, delimiter=" ")
np.savetxt('UMARescueArea_VCM_D2.csv', D2, delimiter=" ")
np.savetxt('UMARescueArea_VCM_AspectX.csv', tri_aspectMapX, delimiter=" ")
np.savetxt('UMARescueArea_VCM_AspectY.csv', tri_aspectMapY, delimiter=" ")

#cmap = plt.cm.gist_earth
#
fig, axes = plt.subplots(constrained_layout=True)
cc = axes.contourf(triX, triY, tri_proximityMap, 100, vmax = 5)
fig.colorbar(cc,location='bottom')
axes.set_aspect('equal')
plt.show()
#
#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(triX, triY, sqrtQ1, 100, cmap = 'Reds')
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()
#
#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(triX, triY, sqrtQ2, 100)
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()
#
#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(triX, triY, D1, 100)
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()
#
#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(triX, triY, D2, 100)
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()

fig, axes = plt.subplots(constrained_layout=True)
cc = axes.contourf(triX, triY, tri_slopeMap, 100, vmax = 30, cmap = 'plasma')
fig.colorbar(cc,location='bottom')
axes.set_aspect('equal')
plt.show()

fig, axes = plt.subplots(constrained_layout=True)
cc = axes.contourf(np.abs(np.arctan2(aspectMapY,aspectMapX)), 20,cmap = 'hsv')
fig.colorbar(cc,location='bottom')
axes.set_aspect('equal')
plt.show()


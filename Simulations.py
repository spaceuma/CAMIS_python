# -*- coding: utf-8 -*-
"""
Simulation Experiment
"""

import numpy as np
import lib.hexgrid as hg
import lib.anisotropic_path_planning as ap
import matplotlib.pyplot as plt
import lib.terrain_processing as tp
import lib.camis as camis
import scipy.interpolate as interp
from time import time
from scipy import ndimage


x = np.linspace(0,50,51)
y = np.linspace(0,50,51)

XX,YY = np.meshgrid(x,y)

#ZZ = np.sin(XX/5)**2 + np.sin(YY/5)**2
#ZZ = 2*np.sin(np.sqrt(XX**2+YY**2)/4)**2
#ZZ = np.sqrt(XX**2+YY**2)/2
#ZZ = XX/5

ZZ = YY/5
ZZ[np.where(YY>40)] = 40/5
ZZ[np.where(YY<20)] = 20/5

slopeMap, aspectMapX, aspectMapY, _ = tp.getConvSlopeMaps(XX, YY, ZZ,\
                                                       1, 1,1)
triX, triY, XY2I, XY2J = hg.getHexGrid(XX,YY,1)
CdRoots,CaRoots,Cl1Roots,Cl2Roots,AniCoLUT = \
camis.readCamis('cuadriga_camis.csv')

CdRoots = [0.0, 0.1, 20.198279377002958]
CaRoots = [0.04681520102276125, 0.8, 20.198279377002958]
Cl1Roots = [2, 20.198279377002958]
Cl2Roots = [2, 20.198279377002958]

maxSlope = np.max(slopeMap)
linearGradient = np.linspace(0,np.ceil(maxSlope),np.ceil(maxSlope)+1)
Cs = []
heading = np.arange(0, 2*np.pi, 0.01)
AniCoLUT = np.zeros((2,linearGradient.size))
Anisotropy = np.zeros_like(linearGradient)
for i,g in enumerate(linearGradient):
    for theta in heading:
        B = camis.computeBeta([1,0] ,theta)
        preCost = camis.computeCAMIScost(B,camis.dirCost(linearGradient[i], CdRoots),\
                                   camis.dirCost(linearGradient[i], CaRoots),\
                                   camis.dirCost(linearGradient[i], Cl1Roots),\
                                   camis.dirCost(linearGradient[i], Cl2Roots))
        Cs.append(preCost)
    Anisotropy[i] = max(Cs)/min(Cs)
    Cs = []
    
AniCoLUT[:][0] = linearGradient
AniCoLUT[:][1] = Anisotropy








obstacleMap = np.zeros_like(slopeMap)
obstacleMap[0,:] = 1
obstacleMap[-1,:] = 1
obstacleMap[:,0] = 1
obstacleMap[:,-1] = 1
proximityMap = ndimage.morphology.distance_transform_edt(1-obstacleMap)
proximityMap[np.where(proximityMap[:]<0)] = 0


init = time()
points = np.zeros((XX.size,2))
points[:,0] = XX.flatten()
points[:,1] = YY.flatten()
slopeMap = interp.griddata(points, slopeMap.flatten(), (triX, triY), method='nearest')
slopeMap[np.where(triX < XX[0,0])] = np.nan
slopeMap[np.where(triX > XX[-1,-1])] = np.nan
slopeMap[np.where(triY < YY[0,0])] = np.nan
slopeMap[np.where(triY > YY[-1,-1])] = np.nan

vectorialData = camis.getVectorialCostMap(\
                                  slopeMap,CdRoots,CaRoots,Cl1Roots,\
                                  Cl2Roots,AniCoLUT)
AspectMap = np.zeros([2,slopeMap.shape[0],slopeMap.shape[1]])
AspectMap[0] = interp.griddata(points, aspectMapX.flatten(), (triX, triY), method='nearest')
AspectMap[1] = interp.griddata(points, aspectMapY.flatten(), (triX, triY), method='nearest')
ProximityMap = interp.griddata(points, proximityMap.flatten(), (triX, triY), method='nearest')
AspectMap[0][np.where(np.isnan(slopeMap))] = np.nan
AspectMap[1][np.where(np.isnan(slopeMap))] = np.nan
ProximityMap[np.where(np.isnan(slopeMap))] = np.nan
obstacleMask = ProximityMap == 0

AnisotropyMap = vectorialData[0][:][:]
AnisotropyMap[obstacleMask] = np.inf

VCMap = np.zeros([4,AnisotropyMap.shape[0],AnisotropyMap.shape[1]])

Q1 = vectorialData[1][:][:]
Q1[obstacleMask] = np.inf
Q2 = vectorialData[2][:][:]
Q2[obstacleMask] = np.inf
D1 = vectorialData[3][:][:]
D1[obstacleMask] = np.inf
D2 = vectorialData[4][:][:]
D2[obstacleMask] = np.inf

VCMap[0] = Q1
VCMap[1] = Q2
VCMap[2] = D1
VCMap[3] = D2

globalRes = triX[0,1]-triX[0,0]
print('Elapsed time to compute the Total Cost Map: '+str(time()-init))



start = np.asarray([40,10])
goal = np.asarray([45,50])

init = time()
Tmap, dirMap, stateMap = \
ap.computeTmap(VCMap,AspectMap,AnisotropyMap,goal,start,triX,triY,globalRes)
print('Elapsed time to compute the Total Cost Map: '+str(time()-init))

IJ2XY = np.zeros([2,triX.shape[0],triY.shape[1]])
IJ2XY[0] = triX
IJ2XY[1] = triY
XY2IJ = np.zeros([2,XY2I.shape[0],XY2J.shape[1]])
XY2IJ[0] = XY2I
XY2IJ[1] = XY2J

startWaypoint = IJ2XY[:,start[1],start[0]]
goalWaypoint = IJ2XY[:,goal[1],goal[0]]
path = ap.getPath(dirMap, IJ2XY, XY2IJ, startWaypoint, goalWaypoint, 0, 0, globalRes)


path = np.asarray(path)

fig, axes = plt.subplots(constrained_layout=True)
cc = axes.contourf(XX,YY,ZZ, 100, cmap = 'plasma')
axes.plot(path[:,0],path[:,1],'r')
fig.colorbar(cc,location='bottom')
axes.set_aspect('equal')
plt.show()


fig, ax = plt.subplots()
ax.contourf(triX, triY, Tmap, 100, cmap = 'magma', alpha = .5)
ax.contour(triX, triY, Tmap, 100, cmap = 'magma')
ax.plot(path[:,0],path[:,1],'r')
ax.set_aspect('equal')
plt.show()

fig, axes = plt.subplots(constrained_layout=True)
cc = axes.contourf(triX, triY, Q1, 100, cmap = 'plasma')
fig.colorbar(cc,location='bottom')
axes.set_aspect('equal')
plt.show()


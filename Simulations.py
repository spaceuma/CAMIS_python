# -*- coding: utf-8 -*-
"""
Simulation Experiment
"""

import numpy as np
import lib.hexgrid as hg
import lib.anisotropic_path_planning as ap
import matplotlib.pyplot as plt

x = np.linspace(0,50,51)
y = np.linspace(0,50,51)



XX,YY = np.meshgrid(x,y)

triX, triY = hg.getHexGrid(XX,YY,1)

globalRes = triX[0,1]-triX[0,0]

AnisotropyMap = np.ones_like(triX)*3
AnisotropyMap[0,:] = np.inf
AnisotropyMap[-1,:] = np.inf
AnisotropyMap[:,0] = np.inf
AnisotropyMap[:,-1] = np.inf

AspectMap = np.zeros([2,AnisotropyMap.shape[0],AnisotropyMap.shape[1]])

aspect = np.pi/4

AspectMap[0], AspectMap[1] = np.ones_like(triX)*np.cos(aspect),\
np.ones_like(triX)*np.sin(aspect)
 
VCMap = np.zeros([4,AnisotropyMap.shape[0],AnisotropyMap.shape[1]])

VCMap[0] = np.ones_like(triX)*2
VCMap[1] = np.ones_like(triX)
VCMap[2] = np.ones_like(triX)
VCMap[3] = np.zeros_like(triX)



goal = np.asarray([20,30])
start = np.asarray([50,30])

TmapG,TmapS, dirMapG, dirMapS, nodeLink, stateMapG,stateMapS = \
ap.computeTmap(VCMap,AspectMap,AnisotropyMap,goal,start,triX,triY,globalRes)

Tmap = np.zeros_like(TmapG)
Tmap[:] = TmapS
for j in range(Tmap.shape[0]):
    for i in range(Tmap.shape[1]):
        if np.isinf(TmapG[j,i]):
            pass
        else:
            Tmap[j,i] = TmapS[nodeLink[1],nodeLink[0]] + TmapG[nodeLink[1],nodeLink[0]] - TmapG[j,i]

fig, ax = plt.subplots()
ax.contourf(triX, triY, Tmap, 20, cmap = 'magma', alpha = .5)
ax.contour(triX, triY, Tmap, 20, cmap = 'magma')
ax.set_aspect('equal')
plt.show()

#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(triX, triY,maxAnisoMap, 100, cmap = 'plasma')
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()


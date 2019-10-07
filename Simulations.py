# -*- coding: utf-8 -*-
"""
Simulation Experiment
"""

import numpy as np
import matplotlib.pyplot as plt
import lib.terrain_processing as tp
import lib.camislib as camislib
import lib.cost_mapping as costmapping
import copy
from time import time


x = np.linspace(0,50,51)
y = np.linspace(0,50,51)

XX,YY = np.meshgrid(x,y)

DEM = 4*(np.cos(XX/4)**2 + np.cos(YY/4)**2)
#DEM = 2*np.sin(np.sqrt((XX/5)**2+(YY/5)**2)/2)**2
#DEM = 2*(np.sin(XX/4)+1.0)
#ZZ = np.sqrt(XX**2+YY**2)/2
#ZZ = XX/5

#DEM = YY/5
#DEM[np.where(YY>35)] = 35/5
#DEM[np.where(YY<15)] = 15/5
demRes = 1.0
planRes = 1.0


# We create a new environment
env1 = costmapping.AnisotropicMap(DEM, demRes, planRes, (0,0))
env1.smoothMap(2.0)
#env1.show3dDEM()

#fig, axes = plt.subplots(constrained_layout=True)
#env1.showMap('slope-deg',fig,axes)
#env1.showMap('laplacian-abs',fig,axes)

env2 = copy.deepcopy(env1)
env3 = copy.deepcopy(env1)
env4 = copy.deepcopy(env1)
#r1 = camislib.CamisModel.fromFile('cuadriga_camis.csv')

cdRoots =  [0.0, 0.20, 1.0]
caRoots =  [0.0, 0.20, 1.0]
cl1Roots = [0.0, 0.20, 1.0]
cl2Roots = [0.0, 0.20, 1.0]

r1 = camislib.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 45.0)
env1.computeVecCostMap(r1)

cdRoots =  [0.0, 0.10, 1.0]
caRoots =  [0.0, 0.30, 1.0]
cl1Roots = [0.0, 0.20, 1.0]
cl2Roots = [0.0, 0.20, 1.0]

r2 = camislib.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 45.0)
env2.computeVecCostMap(r2)

cdRoots =  [0.0, 0.10, 1.0]
caRoots =  [0.0, 0.10, 1.0]
cl1Roots = [0.0, 0.30, 1.0]
cl2Roots = [0.0, 0.10, 1.0]

r3 = camislib.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 45.0)
env3.computeVecCostMap(r3)

cdRoots =  [0.0, 0.10, 1.0]
caRoots =  [0.0, 0.10, 1.0]
cl1Roots = [0.0, 0.10, 1.0]
cl2Roots = [0.0, 0.30, 1.0]

r4 = camislib.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 45.0)
env4.computeVecCostMap(r4)

#start = np.asarray([40,10])
#goal = np.asarray([45,50])
start = np.asarray([10,10])
goal = np.asarray([40,40])

env2.executePlanning(goal,start)
env1.executePlanning(goal,start)
env3.executePlanning(goal,start)
env4.executePlanning(goal,start)

#env2.showResults()

fig, axes = plt.subplots(constrained_layout=True)
env1.showMap('elevation',fig,axes)
env1.showPath(fig,axes,'r')
env2.showPath(fig,axes,'m')
env3.showPath(fig,axes,'g')
env4.showPath(fig,axes,'y')
axes.legend(('CAMIS A', 'CAMIS B', 'CAMIS C', 'CAMIS D'))
axes.set_xlabel('X-axis (m)')
axes.set_ylabel('Y-axis (m)')

fig, axes = plt.subplots(constrained_layout=True)
env1.showPathData('elevation',fig,axes,'r')
env2.showPathData('elevation',fig,axes,'m')
env3.showPathData('elevation',fig,axes,'g')
env4.showPathData('elevation',fig,axes,'y')
plt.grid(True)

fig, axes = plt.subplots(constrained_layout=True)
env1.showPathData('slope',fig,axes,'r')
env2.showPathData('slope',fig,axes,'m')
env3.showPathData('slope',fig,axes,'g')
env4.showPathData('slope',fig,axes,'y')
plt.grid(True)

fig, axes = plt.subplots(constrained_layout=True)
env1.showPathData('pitch',fig,axes,'r')
env2.showPathData('pitch',fig,axes,'m')
env3.showPathData('pitch',fig,axes,'g')
env4.showPathData('pitch',fig,axes,'y')
plt.grid(True)

fig, axes = plt.subplots(constrained_layout=True)
env1.showPathData('roll',fig,axes,'r')
env2.showPathData('roll',fig,axes,'m')
env3.showPathData('roll',fig,axes,'g')
env4.showPathData('roll',fig,axes,'y')
plt.grid(True)

fig, axes = plt.subplots(constrained_layout=True)
env1.showPathData('cost',fig,axes,'r')
env2.showPathData('cost',fig,axes,'m')
env3.showPathData('cost',fig,axes,'g')
env4.showPathData('cost',fig,axes,'y')
plt.grid(True)
axes.set_xlabel('Traversed Distance (m)')
axes.set_ylabel('Cost')

#fig, axes = plt.subplots(constrained_layout=True)
#env1.showPathData('full-orientation',fig,axes,'r')
#env2.showPathData('full-orientation',fig,axes,'m')
#env3.showPathData('full-orientation',fig,axes,'g')
#env4.showPathData('full-orientation',fig,axes,'y')
#plt.grid(True)

#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(XX,YY,ZZ, 100, cmap = 'plasma')
#axes.plot(path[:,0],path[:,1],'r')
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()
#
#
#fig, ax = plt.subplots()
#ax.contourf(triX, triY, Tmap, 100, cmap = 'magma', alpha = .5)
#ax.contour(triX, triY, Tmap, 100, cmap = 'magma')
#ax.plot(path[:,0],path[:,1],'r')
#ax.set_aspect('equal')
#plt.show()
#
#fig, axes = plt.subplots(constrained_layout=True)
#cc = axes.contourf(triX, triY, AnisotropyMap, 100, cmap = 'plasma')
#fig.colorbar(cc,location='bottom')
#axes.set_aspect('equal')
#plt.show()


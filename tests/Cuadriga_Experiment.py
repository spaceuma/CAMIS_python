# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 15:23:37 2019

@author: rsanchez
"""

import numpy as np
import matplotlib.pyplot as plt
import lib.terrain_processing as tp
import lib.camislib as camislib
import lib.cost_mapping as costmapping
import copy
from time import time






hiRes_elevationMap = np.loadtxt(open("terrainData/UMARescueArea_10cmDEM.csv",\
                                     "rb"), delimiter=" ", skiprows=0)
hiRes_posX = np.loadtxt(open("terrainData/UMARescueArea_10cmPosX.csv",\
                                     "rb"), delimiter=" ", skiprows=0)
hiRes_posY = np.loadtxt(open("terrainData/UMARescueArea_10cmPosY.csv",\
                                     "rb"), delimiter=" ", skiprows=0)

hiRes = hiRes_posX[0,1] - hiRes_posX[0,0]

env1 = costmapping.AnisotropicMap(hiRes_elevationMap[600:1400,100:500], hiRes, .2, (hiRes_posX[0,0],hiRes_posY[0,0]))
env1.smoothMap(1.0)

fig, axes = plt.subplots(constrained_layout=True)
env1.showMap('elevation',fig,axes)

cdRoots =  [0.0, 2.0, 20.0]
caRoots =  [0.0, 2.0, 20.0]
cl1Roots = [0.0, 4.0, 20.0]
cl2Roots = [0.0, 4.0, 20.0]

#r1 = camislib.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 45.0)
r1 = camislib.CamisModel.fromFile('cuadriga_camis.csv')
env1.computeVecCostMap(r1)


goal = np.asarray([5,70])
start = np.asarray([27,30])

env1.executePlanning(goal,start)

fig, axes = plt.subplots(constrained_layout=True)
env1.showMap('elevation',fig,axes)
env1.showPath(fig,axes,'r')
# -*- coding: utf-8 -*-
"""
Profiling Planner
"""

import numpy as np
import matplotlib.pyplot as plt
import lib.camislib as camislib
import lib.cost_mapping as costmapping
import copy
from time import time



x = np.linspace(0,50,51)
y = np.linspace(0,50,51)

XX,YY = np.meshgrid(x,y)

#DEM = 4*(np.cos(XX/4)**2 + np.cos(YY/4)**2)
#DEM = 2*np.sin(np.sqrt((XX/5)**2+(YY/5)**2)/2)**2
#DEM = 2*(np.sin(XX/4)+1.0)
#ZZ = np.sqrt(XX**2+YY**2)/2
#ZZ = XX/5

DEM = YY/5
#DEM[np.where(YY>35)] = 35/5
#DEM[np.where(YY<15)] = 15/5
demRes = 1.0
planRes = 1.0


# We create a new environment
env1 = costmapping.AnisotropicMap(DEM, demRes, planRes, (0,0))
env1.smoothMap(1.0)

cdRoots =  [0.0, 0.10, 1.0]
caRoots =  [0.0, 0.20, 1.0]
cl1Roots = [0.0, 0.30, 1.0]
cl2Roots = [0.0, 0.30, 1.0]

r1 = camislib.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 45.0)
env1.computeVecCostMap(r1)

start = np.asarray([10,10])
goal = np.asarray([40,40])

env1.executePlanning(goal,start)
env1.showResults()
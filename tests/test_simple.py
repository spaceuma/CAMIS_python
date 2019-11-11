# -*- coding: utf-8 -*-
"""
Simple test on inclined surface
"""

import numpy as np
from context import camis

x = np.linspace(0,50,51)
y = np.linspace(0,50,51)

XX,YY = np.meshgrid(x,y)

DEM = YY/5
demRes = 1.0
planRes = 1.0


# We create a new environment
env1 = camis.AnisotropicMap(DEM, demRes, planRes, (0,0))
env1.smoothMap(1.0)

cdRoots =  [0.0, 0.10, 1.0]
caRoots =  [0.0, 0.30, 1.0]
cl1Roots = [0.0, 0.20, 1.0]
cl2Roots = [0.0, 0.20, 1.0]

r1 = camis.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 45.0)
env1.computeVecCostMap(r1)

start = np.asarray([10,10])
goal = np.asarray([40,40])

env1.executeBiPlanning(goal,start)

#env1.executePlanning(goal,start)
#env1.showResults()
#
#fig, axes = plt.subplots(constrained_layout=True)
#env1.showMap('elevation',fig,axes)
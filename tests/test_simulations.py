# -*- coding: utf-8 -*-
"""
Simulation Experiment
"""

import numpy as np
import matplotlib.pyplot as plt
from context import camis
import copy

# USER - Choose Simulation Map
nSimMap =4

# =============================================================================
## ENVIRONMENTS CREATION ##
# =============================================================================
x = np.linspace(0,49,50)
y = np.linspace(0,49,50)
XX,YY = np.meshgrid(x,y)


if nSimMap == 1:
    DEM = 4*(np.cos(XX/4)**2 + np.cos(YY/4)**2)
if nSimMap == 2:
    DEM = 2*np.sin(np.sqrt((XX/5)**2+(YY/5)**2)/2)**2
if nSimMap == 3:
    DEM = 2*(np.sin(XX/4)+1.0)
if nSimMap == 4:
    DEM = YY/5
    DEM[np.where(YY>35)] = 35/5
    DEM[np.where(YY<15)] = 15/5

# DEM resolution
demRes = 1.0

# Planning resolution
planRes = 1.0

# We create a new environment (offset is zero,zero) and 4 copies
env1 = camis.AnisotropicMap(DEM, demRes, planRes, (0,0))
env1.smoothMap(2.0)
env2 = copy.deepcopy(env1)
env3 = copy.deepcopy(env1)
env4 = copy.deepcopy(env1)


# =============================================================================
## DIFFERENT CUSTOM CAMIS CREATION ##
# =============================================================================
# ROBOT 1
cdRoots =  [0.0, 0.20, 1.0]
caRoots =  [0.0, 0.20, 1.0]
cl1Roots = [0.0, 0.20, 1.0]
cl2Roots = [0.0, 0.20, 1.0]
r1 = camis.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 90.0)
env1.computeVecCostMap(r1)
r1.showCAMIS()

# ROBOT 2
cdRoots =  [0.0, 0.10, 1.0]
caRoots =  [0.0, 0.30, 1.0]
cl1Roots = [0.0, 0.20, 1.0]
cl2Roots = [0.0, 0.20, 1.0]
r2 = camis.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 90.0)
env2.computeVecCostMap(r2)
r2.showCAMIS()

# ROBOT 3
cdRoots =  [0.0, 0.20, 1.0]
caRoots =  [0.0, 0.20, 1.0]
cl1Roots = [0.0, 0.30, 1.0]
cl2Roots = [0.0, 0.10, 1.0]
r3 = camis.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 90.0)
env3.computeVecCostMap(r3)
r3.showCAMIS()

# ROBOT 4
cdRoots =  [0.0, 0.20, 1.0]
caRoots =  [0.0, 0.20, 1.0]
cl1Roots = [0.0, 0.10, 1.0]
cl2Roots = [0.0, 0.30, 1.0]
r4 = camis.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 90.0)
env4.computeVecCostMap(r4)
r4.showCAMIS()


# =============================================================================
## PATH PLANNING
# =============================================================================

start = np.asarray([10,10])
goal = np.asarray([40,40])

#env1.executeBiPlanning(goal,start)
#env2.executeBiPlanning(goal,start)
#env3.executeBiPlanning(goal,start)
#env4.executeBiPlanning(goal,start)

env1.executePlanning(goal,start)
env2.executePlanning(goal,start)
env3.executePlanning(goal,start)
env4.executePlanning(goal,start)


# =============================================================================
## SHOWING RESULTS
# =============================================================================

fig, axes = plt.subplots(constrained_layout=True)
env1.showMap('elevation',fig,axes)
env1.showPath(fig,axes,'r','solid')
env2.showPath(fig,axes,'b','solid')
env3.showPath(fig,axes,'g','solid')
env4.showPath(fig,axes,'y','solid')
axes.legend(('CAMIS A', 'CAMIS B', 'CAMIS C', 'CAMIS D'))
axes.set_xlabel('X-axis (m)')
axes.set_ylabel('Y-axis (m)')

fig, axes = plt.subplots(constrained_layout=True)
env1.showPathData('elevation',fig,axes,'r')
env2.showPathData('elevation',fig,axes,'b')
env3.showPathData('elevation',fig,axes,'g')
env4.showPathData('elevation',fig,axes,'y')
plt.grid(True)

fig, axes = plt.subplots(2,1,constrained_layout=True)
env1.showPathData('full-orientation',fig,axes[0],'r')
env2.showPathData('full-orientation',fig,axes[0],'b')
env3.showPathData('full-orientation',fig,axes[0],'g')
env4.showPathData('full-orientation',fig,axes[0],'y')
axes[0].legend(('A Slope', 'A Pitch', 'A Roll',\
             'B Slope', 'B Pitch', 'B Roll',
             'C Slope', 'C Pitch', 'C Roll',
             'D Slope', 'D Pitch', 'D Roll'),shadow=True,\
              bbox_to_anchor=(1.05, 1), loc='upper left', \
              borderaxespad=0.,fontsize='x-small')
#axes[0].set_xlabel('X-axis (m)')
axes[0].set_ylabel('Angle (degrees)')
axes[0].grid(True)
env1.showPathData('full-heading',fig,axes[1],'r')
env2.showPathData('full-heading',fig,axes[1],'b')
env3.showPathData('full-heading',fig,axes[1],'g')
env4.showPathData('full-heading',fig,axes[1],'y')
axes[1].legend(('A Heading', 'A Aspect', 'A Beta',\
             'B Heading', 'B Aspect', 'B Beta',
             'C Heading', 'C Aspect', 'C Beta',
             'D Heading', 'D Aspect', 'D Beta'),shadow=True,\
             bbox_to_anchor=(1.05, 1), loc='upper left', \
             borderaxespad=0.,fontsize='x-small')
#axes[1].set_xlabel('X-axis (m)')
axes[1].set_ylabel('Angle (degrees)')
axes[0].get_shared_x_axes().join(axes[0], axes[1])
axes[0].set_xticklabels([])
axes[1].set_xlabel('Traversed Distance (m)')
axes[1].grid(True)


fig, axes = plt.subplots(2,1,constrained_layout=True)
env1.showPathData('cost',fig,axes[0],'r')
env2.showPathData('cost',fig,axes[0],'b')
env3.showPathData('cost',fig,axes[0],'g')
env4.showPathData('cost',fig,axes[0],'y')
axes[0].legend(('A Cost', 'B Cost', 'C Cost', 'D Cost'),\
                shadow=True,bbox_to_anchor=(1.05, 1), loc='upper left', \
                borderaxespad=0.,fontsize='x-small')
#axes[0].set_xlabel('X-axis (m)')
axes[0].set_ylabel('Cost (cost units)')
axes[0].grid(True)
env1.showPathData('total-cost',fig,axes[1],'r')
env2.showPathData('total-cost',fig,axes[1],'b')
env3.showPathData('total-cost',fig,axes[1],'g')
env4.showPathData('total-cost',fig,axes[1],'y')
axes[1].legend(('A Total Cost', 'B Total Cost', 'C Total Cost', \
                'D Total Cost'),shadow=True,bbox_to_anchor=(1.05, 1), \
                loc='upper left', borderaxespad=0.,fontsize='x-small')
#axes[1].set_xlabel('X-axis (m)')
axes[1].set_ylabel('Total Cost (cost units x m)')
axes[0].get_shared_x_axes().join(axes[0], axes[1])
axes[0].set_xticklabels([])
axes[1].set_xlabel('Traversed Distance (m)')
axes[1].grid(True)
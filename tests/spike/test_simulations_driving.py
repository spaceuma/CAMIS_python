"""
MIT License
-----------

Copyright (c) 2021 University of Malaga
Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

Authors: J. Ricardo Sánchez Ibáñez, Carlos J. Pérez del Pulgar
Affiliation: Department of Systems Engineering and Automation
Space Robotics Lab (www.uma.es/space-robotics)
"""
# -*- coding: utf-8 -*-
"""
Simulation Experiment
"""

import numpy as np
import matplotlib.pyplot as plt
from context import camis
import copy

# USER - Choose Simulation Map
nSimMap = 1

# =============================================================================
## ENVIRONMENTS CREATION ##
# =============================================================================
x = np.linspace(0,49,50)
y = np.linspace(0,49,50)
XX,YY = np.meshgrid(x,y)

if nSimMap == 0:
    DEM = 2*(np.cos(XX/8)**2 + np.sin(YY/8)**2)
if nSimMap == 1:
    DEM = (np.cos(XX/4)**2 + np.cos(YY/4)**2)
if nSimMap == 2:
    DEM = 2*np.sin(np.sqrt((XX/5)**2+(YY/5)**2)/2)**2
if nSimMap == 3:
    DEM = 2*(np.sin(XX/4)+1.0)
if nSimMap == 4:
    DEM = YY/5
#    DEM[np.where(YY>35)] = 35/5
#    DEM[np.where(YY<15)] = 15/5
if nSimMap == 5:
    DEM = YY/5 + (np.cos(XX/4)**2 + np.sin(YY/4)**2)

if nSimMap == 6:
    DEM = 2*(np.sin(XX/10+np.pi/4)**2 + np.sin(YY/10+np.pi/4)**2)



# DEM resolution
demRes = 1.0

# Planning resolution
planRes = 1.0

# We create a new environment (offset is zero,zero) and 4 copies
env1 = camis.AnisotropicMap(DEM, demRes, planRes, (0,0))
env1.smoothMap(2.0)
sdThreshold = 11.0
slopeThreshold = 25.0
env1.computeObstacles(sdThreshold,slopeThreshold)
env2 = copy.deepcopy(env1)
env3 = copy.deepcopy(env1)
env4 = copy.deepcopy(env1)

env1.show3dDEM()

# =============================================================================
## DIFFERENT CUSTOM CAMIS CREATION ##
# =============================================================================
# ROBOT 1
r1 = camis.CamisDrivingModel(25.0, .6, .6,  9.8 * 100)
env1.computeVecCostMap(r1)
r1.showCAMIS()

# ROBOT 2
r2 = camis.CamisDrivingModel(25.0, .7, .7,  9.8 * 100)
env2.computeVecCostMap(r2)
r2.showCAMIS()

# ROBOT 3
r3 = camis.CamisDrivingModel(25.0, .8, .8,  9.8 * 100)
env3.computeVecCostMap(r3)
r3.showCAMIS()

# ROBOT 4
r4 = camis.CamisDrivingModel(25.0, .6, .9,  9.8 * 100)
env4.computeVecCostMap(r4)
r4.showCAMIS()


# =============================================================================
## PATH PLANNING
# =============================================================================

goal = np.asarray([10,10])
start = np.asarray([40,40])

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
env1.showMap('proximity',fig,axes)

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
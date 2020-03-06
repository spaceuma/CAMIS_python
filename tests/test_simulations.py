# -*- coding: utf-8 -*-
"""
Simulation Experiment
"""

import numpy as np
import matplotlib.pyplot as plt
from context import camis
import copy
import yaml
import matplotlib.gridspec as gridspec
try:
    from scipy import signal
except:
    raise ImportError('ERROR: scipy module could not be imported')

# USER - Choose Simulation Map
nSimMap = 0

# =============================================================================
## ENVIRONMENTS CREATION ##
# =============================================================================

if nSimMap == 0:
#    x = np.linspace(0,79,80)
#    y = np.linspace(0,79,80)
#    XX,YY = np.meshgrid(x,y)
#    DEM = 0*YY
#    DEM[np.where(YY>15)] = (YY[np.where(YY>15)] - 15)
#    DEM[np.where(YY>25)] = 10
#    DEM[np.where(YY>45)] = 10 + (45-YY[np.where(YY>45)])
#    DEM[np.where(YY>65)] = -10
 
    x = np.linspace(0,49,50)
    y = np.linspace(0,99,100)
    XX,YY = np.meshgrid(x,y)
    DEM = XX
    DEM[np.where(XX>35)] = 35
    DEM[np.where(XX<15)] = 15
    r = 2
    r = r + 1 - r%2
    y,x = np.ogrid[-r: r+1, -r: r+1]
    convMatrix = x**2+y**2 <= r**2
    convMatrix = convMatrix.astype(float)
    DEM = 0.5*signal.convolve2d(DEM, convMatrix/convMatrix.sum(), \
                                      mode='same', boundary='symm') 
    goal = np.asarray([45,95])
    start = np.asarray([5,5])
    
if nSimMap == 1:
    x = np.linspace(0,49,50)
    y = np.linspace(0,99,100)
    XX,YY = np.meshgrid(x,y)
    tXX = np.abs(XX)
    tYY = np.abs(YY - 60)
    dist2center = np.sqrt(tXX**2 + tYY**2)
    DEM1 = np.zeros_like(dist2center)
    DEM1[:] = dist2center   
#    DEM1[np.where(dist2center>30)] = 30
    DEM1 = .5*(DEM1 - 30) + 1.25
    DEM0 = np.ones_like(dist2center)*1.25
    DEM1 = np.maximum(DEM0,DEM1)
#    DEM1 = DEM1 - 30
#    DEM1 = .125*DEM1
    DEM2 = (100-YY)*0.125
#    DEM2[:] = .125*YY - 70/8
    DEM = np.minimum(DEM1,DEM2)
#    DEM[np.where(DEM>0)] = 0
    r = 2
    r = r + 1 - r%2
    y,x = np.ogrid[-r: r+1, -r: r+1]
    convMatrix = x**2+y**2 <= r**2
    convMatrix = convMatrix.astype(float)
    DEM = DEM*(1 + (50-XX)/50*0.2)
    DEM = signal.convolve2d(DEM, convMatrix/convMatrix.sum(), \
                                      mode='same', boundary='symm')
    goal = np.asarray([10,50])
    start = np.asarray([5,5])

if nSimMap == 7:
    DEM = 2*(np.cos(XX/4)**2 + np.sin(YY/4)**2)
if nSimMap == 8:
    DEM = 4*(np.cos(XX/4)**2 + np.cos(YY/4)**2)
if nSimMap == 2:
    DEM = 2*np.sin(np.sqrt((XX/5)**2+(YY/5)**2)/2)**2
if nSimMap == 3:
    DEM = 2*(np.sin(XX/4)+1.0)
if nSimMap == 4:
    DEM = YY/5
    DEM[np.where(YY>35)] = 35/5
    DEM[np.where(YY<15)] = 15/5
if nSimMap == 5:
    DEM = YY/5 + (np.cos(XX/4)**2 + np.sin(YY/4)**2)

if nSimMap == 6:
    x = np.linspace(0,49,50)
    y = np.linspace(0,99,100)
    XX,YY = np.meshgrid(x,y)
    DEM = 5*(np.sin(XX/20)**2 + np.sin(YY/20+np.pi/2)**2)
    


# DEM resolution
demRes = 1.0

# Planning resolution
planRes = 1.0

# We create a new environment (offset is zero,zero) and 4 copies
env1 = camis.AnisotropicMap(DEM, demRes, planRes, (0,0))
#env1.smoothMap(2.0)
#sdThreshold = 11.0
#slopeThreshold = 24.0
#env1.computeObstacles(sdThreshold,slopeThreshold)
env2 = copy.deepcopy(env1)
env3 = copy.deepcopy(env1)
env4 = copy.deepcopy(env1)
env5 = copy.deepcopy(env1)
env6 = copy.deepcopy(env1)

env1.show3dDEM()


fig, axes = plt.subplots(constrained_layout=True)
env1.showMap('slope-deg',fig,axes)
axes.set_xlabel('X-axis (m)')
axes.set_ylabel('Y-axis (m)')

# =============================================================================
## DIFFERENT CUSTOM CAMIS CREATION ##
# =============================================================================

if nSimMap == 0:
    # ROBOT 1
    with open("data/sim02/simRobot1.yml", 'r') as file:
        robot1 = yaml.full_load(file) 
    r1 = camis.CamisDrivingModel(robot1) 
    env1.computeVecCostMap(r1)
    env7 = copy.deepcopy(env1)
    r1.showCAMIS()
    
    # ROBOT 2
    with open("data/sim02/simRobot2.yml", 'r') as file:
        robot2 = yaml.full_load(file)
    r2 = camis.CamisDrivingModel(robot2)
    env2.computeVecCostMap(r2)
    env8 = copy.deepcopy(env2)
    r2.showCAMIS()
    
    # ROBOT 3
    with open("data/sim02/simRobot3.yml", 'r') as file:
        robot3 = yaml.full_load(file)
    r3 = camis.CamisDrivingModel(robot3)
    env3.computeVecCostMap(r3)
    env9 = copy.deepcopy(env3)
    r3.showCAMIS()
    
    # ROBOT 4
    with open("data/sim02/simRobot4.yml", 'r') as file:
        robot4 = yaml.full_load(file)
    r4 = camis.CamisDrivingModel(robot4)
    env4.computeVecCostMap(r4)
    env10 = copy.deepcopy(env4)
    r4.showCAMIS()
    
    # ROBOT 5
    with open("data/sim02/simRobot5.yml", 'r') as file:
        robot5 = yaml.full_load(file)
    r5 = camis.CamisDrivingModel(robot5)
    env5.computeVecCostMap(r5)
    env11 = copy.deepcopy(env5)
    r5.showCAMIS()
    
    # ROBOT 6
    with open("data/sim02/simRobot6.yml", 'r') as file:
        robot6 = yaml.full_load(file)
    r6 = camis.CamisDrivingModel(robot6)
    env6.computeVecCostMap(r6)
    env12 = copy.deepcopy(env6)
    r6.showCAMIS()
    
else:
    
    # ROBOT 1
    with open("data/sim01/simRobot1.yml", 'r') as file:
        robot1 = yaml.full_load(file) 
    r1 = camis.CamisDrivingModel(robot1) 
    env1.computeVecCostMap(r1)
    env7 = copy.deepcopy(env1)
    r1.showCAMIS()
    
    # ROBOT 2
    with open("data/sim01/simRobot2.yml", 'r') as file:
        robot2 = yaml.full_load(file)
    r2 = camis.CamisDrivingModel(robot2)
    env2.computeVecCostMap(r2)
    env8 = copy.deepcopy(env2)
    r2.showCAMIS()
    
    # ROBOT 3
    with open("data/sim01/simRobot3.yml", 'r') as file:
        robot3 = yaml.full_load(file)
    r3 = camis.CamisDrivingModel(robot3)
    env3.computeVecCostMap(r3)
    env9 = copy.deepcopy(env3)
    r3.showCAMIS()
    
    # ROBOT 4
    with open("data/sim01/simRobot4.yml", 'r') as file:
        robot4 = yaml.full_load(file)
    r4 = camis.CamisDrivingModel(robot4)
    env4.computeVecCostMap(r4)
    env10 = copy.deepcopy(env4)
    r4.showCAMIS()
    
    # ROBOT 5
    with open("data/sim01/simRobot5.yml", 'r') as file:
        robot5 = yaml.full_load(file)
    r5 = camis.CamisDrivingModel(robot5)
    env5.computeVecCostMap(r5)
    env11 = copy.deepcopy(env5)
    r5.showCAMIS()
    
    # ROBOT 6
    with open("data/sim01/simRobot6.yml", 'r') as file:
        robot6 = yaml.full_load(file)
    r6 = camis.CamisDrivingModel(robot6)
    env6.computeVecCostMap(r6)
    env12 = copy.deepcopy(env6)
    r6.showCAMIS()

# =============================================================================
## PATH PLANNING
# =============================================================================



env1.executeBiPlanning(goal,start)
env2.executeBiPlanning(goal,start)
env3.executeBiPlanning(goal,start)
env4.executeBiPlanning(goal,start)
env5.executeBiPlanning(goal,start)
env6.executeBiPlanning(goal,start)
env7.executeBiPlanning(start,goal)
env8.executeBiPlanning(start,goal)
env9.executeBiPlanning(start,goal)
env10.executeBiPlanning(start,goal)
env11.executeBiPlanning(start,goal)
env12.executeBiPlanning(start,goal)

#env1.executePlanning(goal,start)
#env2.executePlanning(goal,start)
#env3.executePlanning(goal,start)
#env4.executePlanning(goal,start)
#env5.executePlanning(goal,start)
#env6.executePlanning(goal,start)


# =============================================================================
## SHOWING RESULTS
# =============================================================================
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

fig, axes = plt.subplots(constrained_layout=True)
#fig.suptitle('Simulation 01 - First Run', fontsize=16)
env1.showMap('elevation',fig,axes)
env1.showPath(fig,axes,'r','solid')
env2.showPath(fig,axes,'b','solid')
env3.showPath(fig,axes,'g','solid')
env4.showPath(fig,axes,'y','solid')
env5.showPath(fig,axes,'c','solid')
env6.showPath(fig,axes,'m','solid')
axes.legend(('CAMIS A', 'CAMIS B', 'CAMIS C', 'CAMIS D', 'CAMIS E', 'CAMIS F'))
axes.set_xlabel('X-axis (m)')
axes.set_ylabel('Y-axis (m)')

fig, axes = plt.subplots(constrained_layout=True)
env1.showMap('slope-deg',fig,axes)
env7.showPath(fig,axes,'r','solid')
env8.showPath(fig,axes,'b','solid')
env9.showPath(fig,axes,'g','solid')
env10.showPath(fig,axes,'y','solid')
env11.showPath(fig,axes,'c','solid')
env12.showPath(fig,axes,'m','solid')
axes.legend(('CAMIS A', 'CAMIS B', 'CAMIS C', 'CAMIS D', 'CAMIS E', 'CAMIS F'))
axes.set_xlabel('X-axis (m)')
axes.set_ylabel('Y-axis (m)')



#env1.showPathData('cost',fig,ax2,'r')
#env2.showPathData('cost',fig,ax2,'b')
#env3.showPathData('cost',fig,ax2,'g')
#env4.showPathData('cost',fig,ax2,'y')
#ax2.set_ylabel('C')
#ax2.grid(True)
#
#env1.showPathData('total-cost',fig,ax3,'r')
#env2.showPathData('total-cost',fig,ax3,'b')
#env3.showPathData('total-cost',fig,ax3,'g')
#env4.showPathData('total-cost',fig,ax3,'y')
##axes[2].legend(('A Total Cost', 'B Total Cost', 'C Total Cost', \
##                'D Total Cost'),shadow=True,bbox_to_anchor=(1.05, 1), \
##                loc='upper left', borderaxespad=0.,fontsize='x-small')
##axes[1].set_xlabel('X-axis (m)')
#ax3.set_ylabel('T')
#ax2.get_shared_x_axes().join(ax2, ax3)
#ax2.set_xticklabels([])
#ax3.set_xlabel('Traversed Distance (m)')
#ax3.grid(True)




#fig, axes = plt.subplots(constrained_layout=True)
#env1.showPathData('elevation',fig,axes,'r')
#env2.showPathData('elevation',fig,axes,'b')
#env3.showPathData('elevation',fig,axes,'g')
#env4.showPathData('elevation',fig,axes,'y')
#plt.grid(True)
#
#fig, axes = plt.subplots(2,1,constrained_layout=True)
#env1.showPathData('full-orientation',fig,axes[0],'r')
#env2.showPathData('full-orientation',fig,axes[0],'b')
#env3.showPathData('full-orientation',fig,axes[0],'g')
#env4.showPathData('full-orientation',fig,axes[0],'y')
#axes[0].legend(('A Slope', 'A Pitch', 'A Roll',\
#             'B Slope', 'B Pitch', 'B Roll',
#             'C Slope', 'C Pitch', 'C Roll',
#             'D Slope', 'D Pitch', 'D Roll'),shadow=True,\
#              bbox_to_anchor=(1.05, 1), loc='upper left', \
#              borderaxespad=0.,fontsize='x-small')
##axes[0].set_xlabel('X-axis (m)')
#axes[0].set_ylabel('Angle (degrees)')
#axes[0].grid(True)
#env1.showPathData('full-heading',fig,axes[1],'r')
#env2.showPathData('full-heading',fig,axes[1],'b')
#env3.showPathData('full-heading',fig,axes[1],'g')
#env4.showPathData('full-heading',fig,axes[1],'y')
#axes[1].legend(('A Heading', 'A Aspect', 'A Beta',\
#             'B Heading', 'B Aspect', 'B Beta',
#             'C Heading', 'C Aspect', 'C Beta',
#             'D Heading', 'D Aspect', 'D Beta'),shadow=True,\
#             bbox_to_anchor=(1.05, 1), loc='upper left', \
#             borderaxespad=0.,fontsize='x-small')
##axes[1].set_xlabel('X-axis (m)')
#axes[1].set_ylabel('Angle (degrees)')
#axes[0].get_shared_x_axes().join(axes[0], axes[1])
#axes[0].set_xticklabels([])
#axes[1].set_xlabel('Traversed Distance (m)')
#axes[1].grid(True)
#
#
#fig, axes = plt.subplots(2,1,constrained_layout=True)
#env1.showPathData('cost',fig,axes[0],'r')
#env2.showPathData('cost',fig,axes[0],'b')
#env3.showPathData('cost',fig,axes[0],'g')
#env4.showPathData('cost',fig,axes[0],'y')
#axes[0].legend(('A Cost', 'B Cost', 'C Cost', 'D Cost'),\
#                shadow=True,bbox_to_anchor=(1.05, 1), loc='upper left', \
#                borderaxespad=0.,fontsize='x-small')
##axes[0].set_xlabel('X-axis (m)')
#axes[0].set_ylabel('Cost (cost units)')
#axes[0].grid(True)
#env1.showPathData('total-cost',fig,axes[1],'r')
#env2.showPathData('total-cost',fig,axes[1],'b')
#env3.showPathData('total-cost',fig,axes[1],'g')
#env4.showPathData('total-cost',fig,axes[1],'y')
#axes[1].legend(('A Total Cost', 'B Total Cost', 'C Total Cost', \
#                'D Total Cost'),shadow=True,bbox_to_anchor=(1.05, 1), \
#                loc='upper left', borderaxespad=0.,fontsize='x-small')
##axes[1].set_xlabel('X-axis (m)')
#axes[1].set_ylabel('Total Cost (cost units x m)')
#axes[0].get_shared_x_axes().join(axes[0], axes[1])
#axes[0].set_xticklabels([])
#axes[1].set_xlabel('Traversed Distance (m)')
#axes[1].grid(True)
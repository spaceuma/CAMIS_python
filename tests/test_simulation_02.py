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
nSimMap = 1

# =============================================================================
## ENVIRONMENT CREATION ##
# =============================================================================

x = np.linspace(0,49,50)
y = np.linspace(0,99,100)
XX,YY = np.meshgrid(x,y)
tXX = np.abs(XX)
tYY = np.abs(YY - 60)
dist2center = np.sqrt(tXX**2 + tYY**2)
DEM1 = np.zeros_like(dist2center)
DEM1[:] = dist2center   
DEM1 = .5*(DEM1 - 30) + 1.25
DEM0 = np.ones_like(dist2center)*1.25
DEM1 = np.maximum(DEM0,DEM1)
DEM2 = (100-YY)*0.125
DEM = np.minimum(DEM1,DEM2)
r = 2
r = r + 1 - r%2
y,x = np.ogrid[-r: r+1, -r: r+1]
convMatrix = x**2+y**2 <= r**2
convMatrix = convMatrix.astype(float)
DEM = DEM*(1 + (50-XX)/50*0.2)
DEM = signal.convolve2d(DEM, convMatrix/convMatrix.sum(), \
                                      mode='same', boundary='symm')
goal = np.asarray([10,50])
#start = np.asarray([5,5])
start = np.asarray([45,35])
 


# DEM resolution
demRes = 1.0

# Planning resolution
planRes = 1.0

# We create a new environment (offset is zero,zero) and 4 copies
env1 = camis.AnisotropicMap(DEM, demRes, planRes, (0,0))
env2 = copy.deepcopy(env1)
env3 = copy.deepcopy(env1)
env4 = copy.deepcopy(env1)
env5 = copy.deepcopy(env1)
env6 = copy.deepcopy(env1)

env1.show3dDEM()


fig, axes = plt.subplots(constrained_layout=True)
env1.showMap('elevation',fig,axes)
axes.set_xlabel('X-axis (m)')
axes.set_ylabel('Y-axis (m)')

# =============================================================================
## DIFFERENT CUSTOM CAMIS CREATION ##
# =============================================================================
    
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
plt.style.use('default')
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
env1.showMap('elevation',fig,axes)
env7.showPath(fig,axes,'r','solid')
env8.showPath(fig,axes,'b','solid')
env9.showPath(fig,axes,'g','solid')
env10.showPath(fig,axes,'y','solid')
env11.showPath(fig,axes,'c','solid')
env12.showPath(fig,axes,'m','solid')
axes.legend(('CAMIS A', 'CAMIS B', 'CAMIS C', 'CAMIS D', 'CAMIS E', 'CAMIS F'))
axes.set_xlabel('X-axis (m)')
axes.set_ylabel('Y-axis (m)')

# From https://matplotlib.org/3.2.1/gallery/lines_bars_and_markers/barchart.html#sphx-glr-gallery-lines-bars-and-markers-barchart-py
def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{0:.2f}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
slipCoeffsLabels = ['0.5','1.0','2.0']
integratedTisoGo = [env1.pathComputedTotalCost[-1], env3.pathComputedTotalCost[-1], env5.pathComputedTotalCost[-1]]
integratedTisoReturn = [env7.pathComputedTotalCost[-1], env9.pathComputedTotalCost[-1], env11.pathComputedTotalCost[-1]]
integratedTisoSum = [env1.pathComputedTotalCost[-1] + env7.pathComputedTotalCost[-1],
                     env3.pathComputedTotalCost[-1] + env9.pathComputedTotalCost[-1], 
                     env5.pathComputedTotalCost[-1] + env11.pathComputedTotalCost[-1]]
integratedTanisoGo = [env2.pathComputedTotalCost[-1], env4.pathComputedTotalCost[-1], env6.pathComputedTotalCost[-1]]
integratedTanisoReturn = [env8.pathComputedTotalCost[-1], env10.pathComputedTotalCost[-1], env12.pathComputedTotalCost[-1]]
integratedTanisoSum = [env2.pathComputedTotalCost[-1] + env8.pathComputedTotalCost[-1], 
                       env4.pathComputedTotalCost[-1] + env10.pathComputedTotalCost[-1],
                       env6.pathComputedTotalCost[-1] + env12.pathComputedTotalCost[-1]]
x = np.arange(len(slipCoeffsLabels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, integratedTisoGo, width, label='Isotropic (ρ = 0.8)')
rects1 = ax.bar(x - width/2, integratedTisoReturn, width, bottom = integratedTisoGo, label='Isotropic (ρ = 0.8)')
rects2 = ax.bar(x + width/2, integratedTanisoGo, width, label='Anisotropic (ρ = 0.8)')
rects2 = ax.bar(x + width/2, integratedTanisoReturn, width, bottom = integratedTanisoGo, label='Anisotropic (ρ = 0.8)')
autolabel(rects1)
autolabel(rects2)
ax.set_ylabel('Total Cost [As]')
ax.set_xlabel('Slip Coefficient')
ax.set_xticks(x)
ax.set_xticklabels(slipCoeffsLabels)
ax.legend()
plt.show()


traversedDisoGo = [env1.pathTravDist[-1], env3.pathTravDist[-1], env5.pathTravDist[-1]]
traversedDisoReturn = [env7.pathTravDist[-1], env9.pathTravDist[-1], env11.pathTravDist[-1]]
traversedDanisoGo = [env2.pathTravDist[-1], env4.pathTravDist[-1], env6.pathTravDist[-1]]
traversedDanisoReturn = [env8.pathTravDist[-1], env10.pathTravDist[-1], env12.pathTravDist[-1]]

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, traversedDisoGo, width, label='Isotropic (ρ = 0.8)')
rects1 = ax.bar(x - width/2, traversedDisoReturn, width, bottom = integratedTisoGo, label='Isotropic (ρ = 0.8)')
rects2 = ax.bar(x + width/2, traversedDanisoGo, width, label='Anisotropic (ρ = 0.8)')
rects2 = ax.bar(x + width/2, traversedDanisoReturn, width, bottom = integratedTanisoGo, label='Anisotropic (ρ = 0.8)')
autolabel(rects1)
autolabel(rects2)
ax.set_ylabel('Traversed distance [m]')
ax.set_xlabel('Slip Coefficient')
ax.set_xticks(x)
ax.set_xticklabels(slipCoeffsLabels)
ax.legend()
plt.show()


traversedTimeisoGo = [env1.pathTravDist[-1], env3.pathTravDist[-1], env5.pathTravDist[-1]]
traversedTimeisoReturn = [env7.pathTravDist[-1], env9.pathTravDist[-1], env11.pathTravDist[-1]]
traversedTimeanisoGo = [env2.pathTravDist[-1], env4.pathTravDist[-1], env6.pathTravDist[-1]]
traversedTimeanisoReturn = [env8.pathTravDist[-1], env10.pathTravDist[-1], env12.pathTravDist[-1]]

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, traversedTimeisoGo, width, label='Isotropic (ρ = 0.8)')
rects1 = ax.bar(x - width/2, traversedTimeisoReturn, width, bottom = integratedTisoGo, label='Isotropic (ρ = 0.8)')
rects2 = ax.bar(x + width/2, traversedTimeanisoGo, width, label='Anisotropic (ρ = 0.8)')
rects2 = ax.bar(x + width/2, traversedTimeanisoReturn, width, bottom = integratedTanisoGo, label='Anisotropic (ρ = 0.8)')
autolabel(rects1)
autolabel(rects2)
ax.set_ylabel('Elapsed time [m]')
ax.set_xlabel('Slip Coefficient')
ax.set_xticks(x)
ax.set_xticklabels(slipCoeffsLabels)
ax.legend()
plt.show()


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
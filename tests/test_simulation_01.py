# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 15:23:37 2019

@author: rsanchez
"""

import numpy as np
import matplotlib.pyplot as plt
from context import camis
import copy
import zipfile
import yaml
        
        
# =============================================================================
## SCENE SELECTION
# =============================================================================
# A
# B
# C

# =============================================================================
## LOADING DEM
# =============================================================================
hiRes_elevationMap = np.loadtxt(\
                        open("data/terrainData/UMATerrainCuesta_10cmDEM.csv",\
                        "rb"), delimiter=" ", skiprows=0)
    
hiRes = 0.1
offset = np.loadtxt(open("data/terrainData/UMATerrainCuesta_10cmOffset.csv",\
                                 "rb"), delimiter=" ", skiprows=0)
env = camis.AnisotropicMap(hiRes_elevationMap, hiRes, 0.3,\
                               offset)
posA = np.asarray([5,36]) #Very good
posB = np.asarray([30,60])
posC = np.asarray([30,20])
print('TEST_DEMO: DEM is loaded')


# =============================================================================
## LOADING DIFFERENT CAMIS
# =============================================================================
with open("data/sim01/cuadriga_aniso_01.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_01 = camis.CamisDrivingModel(cuadriga_data)
aniso_01.showDirCosts()
aniso_01.showCAMIS()
with open("data/sim01/cuadriga_iso_01.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
iso_01 = camis.CamisDrivingModel(cuadriga_data)

with open("data/sim01/cuadriga_aniso_02.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_02 = camis.CamisDrivingModel(cuadriga_data)
aniso_02.showDirCosts()
aniso_02.showCAMIS()
with open("data/sim01/cuadriga_iso_02.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
iso_02 = camis.CamisDrivingModel(cuadriga_data)

with open("data/sim01/cuadriga_aniso_03.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_03 = camis.CamisDrivingModel(cuadriga_data)
aniso_03.showDirCosts()
aniso_03.showCAMIS()
with open("data/sim01/cuadriga_iso_03.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
iso_03 = camis.CamisDrivingModel(cuadriga_data)

print('TEST_DEMO: all CAMIS are loaded')



# =============================================================================
## SCENE PROCESSING
# =============================================================================

env_aniso_01_go = copy.deepcopy(env)
env_iso_01_go = copy.deepcopy(env)
env_aniso_02_go = copy.deepcopy(env)
env_iso_02_go = copy.deepcopy(env)
env_aniso_03_go = copy.deepcopy(env)
env_iso_03_go = copy.deepcopy(env)

env_aniso_01_go.computeVecCostMap(aniso_01)
env_aniso_01_back = copy.deepcopy(env_aniso_01_go) 
env_aniso_01_go2 = copy.deepcopy(env_aniso_01_go)
env_aniso_01_back2 = copy.deepcopy(env_aniso_01_go) 
env_iso_01_go.computeVecCostMap(iso_01)
env_iso_01_go2 = copy.deepcopy(env_iso_01_go)

env_aniso_02_go.computeVecCostMap(aniso_02)
env_aniso_02_back = copy.deepcopy(env_aniso_02_go) 
env_aniso_02_go2 = copy.deepcopy(env_aniso_02_go) 
env_aniso_02_back2 = copy.deepcopy(env_aniso_02_go) 
env_iso_02_go.computeVecCostMap(iso_02)
env_iso_02_go2 = copy.deepcopy(env_iso_02_go)

env_aniso_03_go.computeVecCostMap(aniso_03)
env_aniso_03_back = copy.deepcopy(env_aniso_03_go) 
env_aniso_03_go2 = copy.deepcopy(env_aniso_03_go) 
env_aniso_03_back2 = copy.deepcopy(env_aniso_03_go) 
env_iso_03_go.computeVecCostMap(iso_03)
env_iso_03_go2 = copy.deepcopy(env_iso_03_go)
print('TEST_DEMO: the environments are processed')

# =============================================================================
## EXECUTING PATH PLANNING
# =============================================================================

env_aniso_01_go.executeBiPlanning(posB,posA)
env_aniso_01_back.executeBiPlanning(posA,posB)
env_aniso_01_go2.executeBiPlanning(posC,posA)
env_aniso_01_back2.executeBiPlanning(posA,posC)
env_iso_01_go.executeBiPlanning(posB,posA)
env_iso_01_go2.executeBiPlanning(posC,posA)

env_aniso_02_go.executeBiPlanning(posB,posA)
env_aniso_02_back.executeBiPlanning(posA,posB)
env_aniso_02_go2.executeBiPlanning(posC,posA)
env_aniso_02_back2.executeBiPlanning(posA,posC)
env_iso_02_go.executeBiPlanning(posB,posA)
env_iso_02_go2.executeBiPlanning(posC,posA)

env_aniso_03_go.executeBiPlanning(posB,posA)
env_aniso_03_back.executeBiPlanning(posA,posB)
env_aniso_03_go2.executeBiPlanning(posC,posA)
env_aniso_03_back2.executeBiPlanning(posA,posC)
env_iso_03_go.executeBiPlanning(posB,posA)
env_iso_03_go2.executeBiPlanning(posC,posA)

# =============================================================================
## SHOWING RESULTS
# =============================================================================

env.show3dDEM()

env_aniso_01_go.showHexAnisotropyMap()

fig, axes = plt.subplots(constrained_layout=True)
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
env_aniso_01_go.showMap('hex-elevation',fig,axes)

env_aniso_01_go.showPath(fig,axes,'b','solid')
env_aniso_01_back.showPath(fig,axes,'b','dashed')
env_iso_01_go.showPath(fig,axes,'m','solid')

env_aniso_02_go.showPath(fig,axes,'c','solid')
env_aniso_02_back.showPath(fig,axes,'c','dashed')
env_iso_02_go.showPath(fig,axes,'r','solid')

env_aniso_03_go.showPath(fig,axes,'k','solid')
env_aniso_03_back.showPath(fig,axes,'k','dashed')
env_iso_03_go.showPath(fig,axes,'y','solid')

env_aniso_01_go2.showPath(fig,axes,'b','solid')
env_aniso_01_back2.showPath(fig,axes,'b','dashed')
env_iso_01_go2.showPath(fig,axes,'m','solid')

env_aniso_02_go2.showPath(fig,axes,'c','solid')
env_aniso_02_back2.showPath(fig,axes,'c','dashed')
env_iso_02_go2.showPath(fig,axes,'r','solid')

env_aniso_03_go2.showPath(fig,axes,'k','solid')
env_aniso_03_back2.showPath(fig,axes,'k','dashed')
env_iso_03_go2.showPath(fig,axes,'y','solid')

axes.set_xlabel('X-axis [m]')
axes.set_ylabel('Y-axis [m]')

axes.legend(('CAMIS A(Go)', 'CAMIS A (Return)',\
             'Isotropic A (Go and Return)', \
             'CAMIS B(Go)', 'CAMIS B (Return)', \
             'Isotropic B (Go and Return)', \
             'CAMIS C(Go)', 'CAMIS C (Return)', \
             'Isotropic C (Go and Return)'))

#env_max_back.showPath(fig,axes,'g','solid')
#axes.legend(('aniso_imu (Back)', 'iso_average (Back)',\
#             'iso_maximum (Back)'))

#print('TEST_DEMO: total cost expected by aniso_dem is ' + \
#      str(env_dem_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
#      str(env_dem_back.pathEstimatedTotalCost[0]/3600)+ ' Ah')
print('TEST_DEMO: total cost expected by aniso_imu is ' + \
      str(env_aniso_01_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
      str(env_aniso_01_back.pathEstimatedTotalCost[0]/3600)+ ' Ah')
print('TEST_DEMO: total cost expected by iso_nom is ' + \
      str(env_iso_01_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
      str(env_iso_01_go.pathEstimatedTotalCost[0]/3600)+ ' Ah')
print('TEST_DEMO: total cost expected by aniso_imu is ' + \
      str(env_aniso_02_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
      str(env_aniso_02_back.pathEstimatedTotalCost[0]/3600)+ ' Ah')
print('TEST_DEMO: total cost expected by iso_nom is ' + \
      str(env_iso_02_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
      str(env_iso_02_go.pathEstimatedTotalCost[0]/3600)+ ' Ah')
#print('TEST_DEMO: total cost expected by iso_average is ' + \
#      str(env_aver_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
#      str(env_aver_back.pathEstimatedTotalCost[0]/3600)+ ' Ah')
#print('TEST_DEMO: total cost expected by iso_maximum is ' + \
#      str(env_max_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
#      str(env_max_back.pathEstimatedTotalCost[0]/3600)+ ' Ah')
#fig, axes = plt.subplots(constrained_layout=True)
#envGo_dem.showPathData('cost',fig,axes,'r')
#envBack_dem.showPathData('cost',fig,axes,'b')
#envGo_imu.showPathData('cost',fig,axes,'m')
#envBack_imu.showPathData('cost',fig,axes,'g')
#axes.legend(('aniso_dem (Go)', 'aniso_dem (Back)', 'aniso_imu (Go)', 'aniso_imu (Back)'))
#axes.set_xlabel('2D-Path length (m)')
#axes.set_ylabel('Cost (A)')
#
plt.style.use('seaborn-darkgrid')
fig, axes = plt.subplots(constrained_layout=True)
env_aniso_01_go.showPathData('total-cost',fig,axes,'b','solid')
env_aniso_01_back.showPathData('total-cost',fig,axes,'b','dashed')
env_iso_01_go.showPathData('total-cost',fig,axes,'m','solid')
env_aniso_02_go.showPathData('total-cost',fig,axes,'c','solid')
env_aniso_02_back.showPathData('total-cost',fig,axes,'c','dashed')
env_iso_02_go.showPathData('total-cost',fig,axes,'r','solid')
plt.style.use('default')

plt.style.use('seaborn-darkgrid')
fig, axes = plt.subplots(constrained_layout=True)
env_aniso_01_go.showPathData('pitch',fig,axes,'b','solid')
env_aniso_01_back.showPathData('pitch',fig,axes,'b','dashed')
env_iso_01_go.showPathData('pitch',fig,axes,'m','solid')
#env_iso_01_back.showPathData('pitch',fig,axes,'m','dashed')
env_aniso_02_go.showPathData('pitch',fig,axes,'c','solid')
env_aniso_02_back.showPathData('pitch',fig,axes,'c','dashed')
env_iso_02_go.showPathData('pitch',fig,axes,'r','solid')
env_aniso_03_go.showPathData('pitch',fig,axes,'k','solid')
env_aniso_03_back.showPathData('pitch',fig,axes,'k','dashed')
env_iso_03_go.showPathData('pitch',fig,axes,'y','solid')
plt.style.use('default')
axes.set_xlabel('Traversed distance [m]')
axes.set_ylabel('Pitch [degrees]')
#axes.legend(('CAMIS A(Go)', 'CAMIS A (Return)',\
#             'Isotropic A (Go and Return)', \
#             'CAMIS B(Go)', 'CAMIS B (Return)', \
#             'Isotropic B (Go and Return)', \
#             'CAMIS C(Go)', 'CAMIS C (Return)', \
#             'Isotropic C (Go and Return)'))

plt.style.use('seaborn-darkgrid')
fig, axes = plt.subplots(constrained_layout=True)
env_aniso_01_go.showPathData('roll',fig,axes,'b','solid')
env_aniso_01_back.showPathData('roll',fig,axes,'b','dashed')
env_iso_01_go.showPathData('roll',fig,axes,'m','solid')
env_aniso_02_go.showPathData('roll',fig,axes,'c','solid')
env_aniso_02_back.showPathData('roll',fig,axes,'c','dashed')
env_iso_02_go.showPathData('roll',fig,axes,'r','solid')
plt.style.use('default')

plt.style.use('seaborn-darkgrid')
optList = ['segment', 'heading']
for i,opt in enumerate(optList):
    fig, axes = plt.subplots(constrained_layout=True)
    env_aniso_01_go.showPathData(opt,fig,axes,'b','solid')
    env_aniso_01_back.showPathData(opt,fig,axes,'b','dashed')
    env_iso_01_go.showPathData(opt,fig,axes,'m','solid')
    env_aniso_02_go.showPathData(opt,fig,axes,'c','solid')
    env_aniso_02_back.showPathData(opt,fig,axes,'c','dashed')
    env_iso_02_go.showPathData(opt,fig,axes,'r','solid')
plt.style.use('default')

plt.style.use('seaborn-darkgrid')
fig, axes = plt.subplots(constrained_layout=True)
env_aniso_01_go.showPathData('beta',fig,axes,'b','solid')
env_aniso_01_back.showPathData('beta',fig,axes,'b','dashed')
env_iso_01_go.showPathData('beta',fig,axes,'m','solid')
env_aniso_02_go.showPathData('beta',fig,axes,'c','solid')
env_aniso_02_back.showPathData('beta',fig,axes,'c','dashed')
env_iso_02_go.showPathData('beta',fig,axes,'r','solid')
plt.style.use('default')


#envBack_bidem.showPathData('total-cost-estimated',fig,axes,'b')
#envGo_biimu.showPathData('total-cost-estimated',fig,axes,'m')
#envBack_biimu.showPathData('total-cost-estimated',fig,axes,'g')

#envGo.show3dDEM()
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

# =============================================================================
## SCENE SELECTION
# =============================================================================
# A
# B
# C
terrain = 'A'

if (terrain != 'A') and (terrain != 'B') and (terrain != 'C'):
    raise NameError('The terrain selected is not valid')

# =============================================================================
## LOADING DEM
# =============================================================================
zipRef = zipfile.ZipFile("data/terrainData/UMARescueArea_10cmDEM.zip","r")
zipRef.extractall("data/terrainData/")
hiRes_elevationMap = np.loadtxt(\
                        open("data/terrainData/UMARescueArea_10cmDEM.csv",\
                             "rb"), delimiter=" ", skiprows=0)
hiRes_posX = np.loadtxt(open("data/terrainData/UMARescueArea_10cmPosX.csv",\
                                     "rb"), delimiter=" ", skiprows=0)
hiRes_posY = np.loadtxt(open("data/terrainData/UMARescueArea_10cmPosY.csv",\
                                     "rb"), delimiter=" ", skiprows=0)
hiRes = hiRes_posX[0,1] - hiRes_posX[0,0]
print('TEST_DEMO: DEM is loaded')

# =============================================================================
## LOADING DIFFERENT CAMIS
# =============================================================================
aniso_dem = camis.CamisModel.fromFile('cuadriga_camis_dem.csv')
aniso_imu = camis.CamisModel.fromFile('cuadriga_camis_imu.csv') 
iso_average = camis.CamisModel.fromFile('cuadriga_camis_dem_iso_med.csv')
iso_maximum = camis.CamisModel.fromFile('cuadriga_camis_dem_iso_asc.csv')
print('TEST_DEMO: all CAMIS are loaded')

# =============================================================================
## CREATING THE ENVIRONMENT
# =============================================================================

if terrain == 'A':
    env = camis.AnisotropicMap(hiRes_elevationMap[600:1100,0:600], hiRes, 0.5,\
                               (hiRes_posX[0,0],hiRes_posY[0,0]))
    posA = np.asarray([10,40])
    posB = np.asarray([50,45])
if terrain == 'B':
    env = camis.AnisotropicMap(hiRes_elevationMap[450:1250,700:1100], hiRes,\
                               0.5, (hiRes_posX[0,0],hiRes_posY[0,0]))
    posA = np.asarray([5,50])
    posB = np.asarray([30,40])
if terrain == 'C':
    env = camis.AnisotropicMap(hiRes_elevationMap[0:550,800:1200], hiRes, 0.5,\
                               (hiRes_posX[0,0],hiRes_posY[0,0]))
    posA = np.asarray([30,50])
    posB = np.asarray([15,5])
print('TEST_DEMO: the environment is set up')

# =============================================================================
## SCENE PROCESSING
# =============================================================================

sdThreshold = 10.0 # Standard Deviation Threshold for roughness avoidance
env.smoothMap(2.0)
env.processSlope(sdThreshold)

env_dem_go = copy.deepcopy(env)
env_imu_go = copy.deepcopy(env)
env_aver_go = copy.deepcopy(env)
env_max_go = copy.deepcopy(env)

env_dem_go.computeVecCostMap(aniso_dem)
env_dem_back = copy.deepcopy(env_dem_go)
env_imu_go.computeVecCostMap(aniso_imu)
env_imu_back = copy.deepcopy(env_imu_go) 
env_aver_go.computeVecCostMap(iso_average)
env_aver_back = copy.deepcopy(env_aver_go)
env_max_go.computeVecCostMap(iso_maximum)
env_max_back = copy.deepcopy(env_max_go)
print('TEST_DEMO: the environment is processed')

# =============================================================================
## EXECUTING PATH PLANNING
# =============================================================================

env_dem_go.executeBiPlanning(posB,posA)
env_dem_back.executeBiPlanning(posA,posB)
env_imu_go.executeBiPlanning(posB,posA)
env_imu_back.executeBiPlanning(posA,posB)
env_aver_go.executeBiPlanning(posB,posA)
env_aver_back.executeBiPlanning(posA,posB)
env_max_go.executeBiPlanning(posB,posA)
env_max_back.executeBiPlanning(posA,posB)


# =============================================================================
## SHOWING RESULTS
# =============================================================================

env.show3dDEM()

fig, axes = plt.subplots(constrained_layout=True)
env.showMap('old-elevation',fig,axes)
env_dem_go.showPath(fig,axes,'r','solid')
env_dem_back.showPath(fig,axes,'r','dotted')
env_imu_go.showPath(fig,axes,'b','solid')
env_imu_back.showPath(fig,axes,'b','dotted')
env_aver_go.showPath(fig,axes,'m','solid')
env_aver_back.showPath(fig,axes,'m','dotted')
env_max_go.showPath(fig,axes,'g','solid')
env_max_back.showPath(fig,axes,'g','dotted')
axes.legend(('aniso_dem (Go)', 'aniso_dem (Back)', 'aniso_imu (Go)', \
             'aniso_imu (Back)', 'iso_average (Go)', 'iso_average (Back)',\
             'iso_maximum (Go)', 'iso_maximum (Back)'))

#fig, axes = plt.subplots(constrained_layout=True)
#envGo_dem.showPathData('cost',fig,axes,'r')
#envBack_dem.showPathData('cost',fig,axes,'b')
#envGo_imu.showPathData('cost',fig,axes,'m')
#envBack_imu.showPathData('cost',fig,axes,'g')
#axes.legend(('aniso_dem (Go)', 'aniso_dem (Back)', 'aniso_imu (Go)', 'aniso_imu (Back)'))
#axes.set_xlabel('2D-Path length (m)')
#axes.set_ylabel('Cost (A)')
#
#fig, axes = plt.subplots(constrained_layout=True)
#envGo_bidem.showPathData('total-cost-estimated',fig,axes,'r')
#envBack_bidem.showPathData('total-cost-estimated',fig,axes,'b')
#envGo_biimu.showPathData('total-cost-estimated',fig,axes,'m')
#envBack_biimu.showPathData('total-cost-estimated',fig,axes,'g')

#envGo.show3dDEM()
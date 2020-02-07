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
try:
    from scipy import signal
except:
    raise ImportError('ERROR: scipy module could not be imported')

# =============================================================================
## SCENE SELECTION
# =============================================================================
# A
# B
# C
terrain = 'C'

if (terrain != 'A') and (terrain != 'B') and (terrain != 'C'):
    raise NameError('The terrain selected is not valid')

# =============================================================================
## LOADING DEM
# =============================================================================
zipRef = zipfile.ZipFile("data/terrainData/UMARescueArea_10cmDEM.zip","r")
zipRef.extractall("data/terrainData/")
hiRes_elevationMap = np.loadtxt(\
                        open("data/terrainData/UMATerrain_10cmDEM.csv",\
                             "rb"), delimiter=" ", skiprows=0)
#hiRes_posX = np.loadtxt(open("data/terrainData/UMARescueArea_10cmPosX.csv",\
#                                     "rb"), delimiter=" ", skiprows=0)
#hiRes_posY = np.loadtxt(open("data/terrainData/UMARescueArea_10cmPosY.csv",\
#                                     "rb"), delimiter=" ", skiprows=0)
hiRes = 0.1
offset = np.loadtxt(\
                        open("data/terrainData/UMATerrain_10cmOffset.csv",\
                             "rb"), delimiter=" ", skiprows=0)
print('TEST_DEMO: DEM is loaded')

# =============================================================================
## LOADING DIFFERENT CAMIS
# =============================================================================

#cdRoots =  [0.0, 0.1, 1.0]
#caRoots =  [0.0, 0.3, 1.0]
#cl1Roots = [0.0, 0.8, 1.0]
#cl2Roots = [0.0, 0.8, 1.0]
#aniso_imu = camis.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 25.0)
#
#cdRoots =  [0.0, 0.40, 1.0]
#caRoots =  [0.0, 0.40, 1.0]
#cl1Roots = [0.0, 0.40, 1.0]
#cl2Roots = [0.0, 0.40, 1.0]
#iso_average = camis.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 25.0)
#
#cdRoots =  [0.0, 0.08, 1.0]
#caRoots =  [0.0, 0.08, 1.0]
#cl1Roots = [0.0, 0.08, 1.0]
#cl2Roots = [0.0, 0.08, 1.0]
#iso_maximum = camis.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 25.0)

#aniso_imu = camis.CamisDrivingModel(25.0, 13.04883169, 14.36578031,  9.6977448)

#aniso_imu = camis.CamisDrivingModel(25.0, .7, .7,  9.8 * 50)
with open("data/cuadriga.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_imu = camis.CamisDrivingModel(cuadriga_data)

with open("data/cuadriga_iso_nominal.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
iso_nominal = camis.CamisDrivingModel(cuadriga_data)

#aniso_imu = camis.CamisModel.fromFile('data/camisRoots/cuadriga_camis_imu.csv') 
#iso_average = camis.CamisModel.fromFile('data/camisRoots/cuadriga_camis_imu_iso_med.csv')
#iso_maximum = camis.CamisModel.fromFile('data/camisRoots/cuadriga_camis_imu_iso_asc.csv')
print('TEST_DEMO: all CAMIS are loaded')

# =============================================================================
## CREATING THE ENVIRONMENT
# =============================================================================

if terrain == 'A':
    newMap = np.zeros_like(hiRes_elevationMap)
    newMap[:] = hiRes_elevationMap
    newMap[1600:1800,700:850] = 60.0
    env = camis.AnisotropicMap(newMap[1000:2600,400:1200], hiRes, 0.5,\
                               offset)
    posA = np.asarray([30,20])
    posB = np.asarray([50,150])
if terrain == 'B':
    DEM = hiRes_elevationMap[1200:1600,550:950]
#    r = 4
#    r = r + 1 - r%2
#    y,x = np.ogrid[-r: r+1, -r: r+1]
#    convMatrix = x**2+y**2 <= r**2
#    convMatrix = convMatrix.astype(float)
#    DEM = 0.5*signal.convolve2d(DEM, convMatrix/convMatrix.sum(), \
#                                      mode='same', boundary='symm')
    env = camis.AnisotropicMap(DEM, hiRes, 0.4,\
                               offset)
    posA = np.asarray([25,10])
    posB = np.asarray([25,35])
if terrain == 'C':
    env = camis.AnisotropicMap(hiRes_elevationMap[800:1400,1100:1500], hiRes, 0.4,\
                               offset)
#    posA = np.asarray([5,60])
#    posB = np.asarray([30,30])
    posA = np.asarray([10,40])
    posB = np.asarray([30,20])
print('TEST_DEMO: the environment is set up')

# =============================================================================
## SCENE PROCESSING
# =============================================================================


env_dem_go = copy.deepcopy(env)
env_imu_go = copy.deepcopy(env)
env_nom_go = copy.deepcopy(env)
env_max_go = copy.deepcopy(env)

env_imu_go.computeVecCostMap(aniso_imu)
env_imu_back = copy.deepcopy(env_imu_go) 
env_nom_go.computeVecCostMap(iso_nominal)
env_nom_back = copy.deepcopy(env_nom_go)

#env_max_go.computeVecCostMap(iso_maximum)
#env_max_back = copy.deepcopy(env_max_go)
print('TEST_DEMO: the environment is processed')

# =============================================================================
## EXECUTING PATH PLANNING
# =============================================================================

#env_dem_go.executePlanning(posB,posA)
#env_dem_back.executePlanning(posA,posB)
env_imu_go.executePlanning(posB,posA)
env_imu_back.executePlanning(posA,posB)
env_nom_go.executePlanning(posB,posA)
env_nom_back.executePlanning(posA,posB)
#env_aver_go.executePlanning(posB,posA)
#env_aver_back.executePlanning(posA,posB)
#env_max_go.executePlanning(posB,posA)
#env_max_back.executePlanning(posA,posB)


# =============================================================================
## SHOWING RESULTS
# =============================================================================

env.show3dDEM()

fig, axes = plt.subplots(constrained_layout=True)
env_imu_go.showMap('standard-deviation',fig,axes)

fig, axes = plt.subplots(constrained_layout=True)
env_imu_go.showMap('slope-deg',fig,axes)

fig, axes = plt.subplots(constrained_layout=True)
env_imu_go.showMap('proximity',fig,axes)

env_imu_go.showHexAnisotropyMap()

fig, axes = plt.subplots(constrained_layout=True)
env.showMap('elevation',fig,axes)
env_imu_go.showPath(fig,axes,'b','solid')
env_nom_go.showPath(fig,axes,'m','solid')
#env_max_go.showPath(fig,axes,'g','solid')
#axes.legend(('aniso_imu (Go)', 'iso_average (Go)', \
#             'iso_maximum (Go)'))
fig, axes = plt.subplots(constrained_layout=True)
env.showMap('elevation',fig,axes)
env_imu_back.showPath(fig,axes,'b','solid')
env_nom_back.showPath(fig,axes,'m','solid')
#env_max_back.showPath(fig,axes,'g','solid')
#axes.legend(('aniso_imu (Back)', 'iso_average (Back)',\
#             'iso_maximum (Back)'))

fig, axes = plt.subplots(constrained_layout=True)
env.showMap('slope-deg',fig,axes)

#print('TEST_DEMO: total cost expected by aniso_dem is ' + \
#      str(env_dem_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
#      str(env_dem_back.pathEstimatedTotalCost[0]/3600)+ ' Ah')
print('TEST_DEMO: total cost expected by aniso_imu is ' + \
      str(env_imu_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
      str(env_imu_back.pathEstimatedTotalCost[0]/3600)+ ' Ah')
print('TEST_DEMO: total cost expected by iso_nom is ' + \
      str(env_nom_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
      str(env_nom_back.pathEstimatedTotalCost[0]/3600)+ ' Ah')
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
#fig, axes = plt.subplots(constrained_layout=True)
#envGo_bidem.showPathData('total-cost-estimated',fig,axes,'r')
#envBack_bidem.showPathData('total-cost-estimated',fig,axes,'b')
#envGo_biimu.showPathData('total-cost-estimated',fig,axes,'m')
#envBack_biimu.showPathData('total-cost-estimated',fig,axes,'g')

#envGo.show3dDEM()
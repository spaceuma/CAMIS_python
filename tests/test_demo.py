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

with open("data/cuadriga_aniso_norisk.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_norisk = camis.CamisDrivingModel(cuadriga_data)
aniso_norisk.showDirCosts()
aniso_norisk.showCAMIS()

with open("data/cuadriga_iso_norisk.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
iso_norisk = camis.CamisDrivingModel(cuadriga_data)

with open("data/cuadriga_aniso_wrisk.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_wrisk = camis.CamisDrivingModel(cuadriga_data)

with open("data/cuadriga_iso_wrisk.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
iso_wrisk = camis.CamisDrivingModel(cuadriga_data)

print('TEST_DEMO: all CAMIS are loaded')


# =============================================================================
## CREATING THE ENVIRONMENT
# =============================================================================

if terrain == 'A':
    env = camis.AnisotropicMap(hiRes_elevationMap[1450:1750,1150:1450], hiRes, 0.4,\
                               offset)
#    posA = np.asarray([5,25]) #Very good
#    posB = np.asarray([25,5])
    
    posA = np.asarray([5,5])
    posB = np.asarray([20,5])
    
#    posA = np.asarray([10,40])
#    posB = np.asarray([30,20])
    
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
    env = camis.AnisotropicMap(hiRes_elevationMap[850:1250,1150:1450], hiRes, 0.4,\
                               offset)
    posA = np.asarray([5,25]) #Very good
#    posB = np.asarray([25,5])
    
#    posA = np.asarray([5,10])
    posB = np.asarray([25,25])
    
#    posA = np.asarray([10,40])
#    posB = np.asarray([30,20])
print('TEST_DEMO: the environment is set up')

# =============================================================================
## SCENE PROCESSING
# =============================================================================

env_aniso_norisk_go = copy.deepcopy(env)
env_iso_norisk_go = copy.deepcopy(env)
env_aniso_wrisk_go = copy.deepcopy(env)
env_iso_wrisk_go = copy.deepcopy(env)

env_aniso_norisk_go.computeVecCostMap(aniso_norisk)
env_aniso_norisk_back = copy.deepcopy(env_aniso_norisk_go) 
env_iso_norisk_go.computeVecCostMap(iso_norisk)
env_iso_norisk_back = copy.deepcopy(env_iso_norisk_go)
env_aniso_wrisk_go.computeVecCostMap(aniso_wrisk)
env_aniso_wrisk_back = copy.deepcopy(env_aniso_wrisk_go) 
env_iso_wrisk_go.computeVecCostMap(iso_wrisk)
env_iso_wrisk_back = copy.deepcopy(env_iso_wrisk_go)

print('TEST_DEMO: the environments are processed')

# =============================================================================
## EXECUTING PATH PLANNING
# =============================================================================

env_aniso_norisk_go.executeBiPlanning(posB,posA)
env_aniso_norisk_back.executeBiPlanning(posA,posB)
env_iso_norisk_go.executeBiPlanning(posB,posA)
env_iso_norisk_back.executeBiPlanning(posA,posB)
env_aniso_wrisk_go.executeBiPlanning(posB,posA)
env_aniso_wrisk_back.executeBiPlanning(posA,posB)
env_iso_wrisk_go.executeBiPlanning(posB,posA)
env_iso_wrisk_back.executeBiPlanning(posA,posB)

# =============================================================================
## SHOWING RESULTS
# =============================================================================

env.show3dDEM()

fig, axes = plt.subplots(constrained_layout=True)
env_aniso_norisk_go.showMap('standard-deviation',fig,axes)

fig, axes = plt.subplots(constrained_layout=True)
env_aniso_norisk_go.showMap('aspect-deg',fig,axes)

fig, axes = plt.subplots(constrained_layout=True)
env_iso_norisk_go.showMap('slope-deg',fig,axes)

env_aniso_norisk_go.showHexAnisotropyMap()

fig, axes = plt.subplots(constrained_layout=True)
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
env.showMap('elevation',fig,axes)
env_aniso_norisk_go.showPath(fig,axes,'b','solid')
env_aniso_norisk_back.showPath(fig,axes,'b','dashed')
env_iso_norisk_go.showPath(fig,axes,'m','solid')
env_iso_norisk_back.showPath(fig,axes,'m','dashed')
env_aniso_wrisk_go.showPath(fig,axes,'c','solid')
env_aniso_wrisk_back.showPath(fig,axes,'c','dashed')
env_iso_wrisk_go.showPath(fig,axes,'r','solid')
env_iso_wrisk_back.showPath(fig,axes,'r','dashed')
axes.set_xlabel('X-axis [m]')
axes.set_ylabel('Y-axis [m]')

#env_max_go.showPath(fig,axes,'g','solid')
#axes.legend(('aniso_imu (Go)', 'iso_average (Go)', \
#             'iso_maximum (Go)'))

#env_max_back.showPath(fig,axes,'g','solid')
#axes.legend(('aniso_imu (Back)', 'iso_average (Back)',\
#             'iso_maximum (Back)'))

fig, axes = plt.subplots(constrained_layout=True)
env.showMap('slope-deg',fig,axes)

#print('TEST_DEMO: total cost expected by aniso_dem is ' + \
#      str(env_dem_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
#      str(env_dem_back.pathEstimatedTotalCost[0]/3600)+ ' Ah')
print('TEST_DEMO: total cost expected by aniso_imu is ' + \
      str(env_aniso_norisk_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
      str(env_aniso_norisk_back.pathEstimatedTotalCost[0]/3600)+ ' Ah')
print('TEST_DEMO: total cost expected by iso_nom is ' + \
      str(env_iso_norisk_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
      str(env_iso_norisk_back.pathEstimatedTotalCost[0]/3600)+ ' Ah')
print('TEST_DEMO: total cost expected by aniso_imu is ' + \
      str(env_aniso_wrisk_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
      str(env_aniso_wrisk_back.pathEstimatedTotalCost[0]/3600)+ ' Ah')
print('TEST_DEMO: total cost expected by iso_nom is ' + \
      str(env_iso_wrisk_go.pathEstimatedTotalCost[0]/3600) + ' Ah and ' + \
      str(env_iso_wrisk_back.pathEstimatedTotalCost[0]/3600)+ ' Ah')
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
env_aniso_norisk_go.showPathData('total-cost',fig,axes,'b','solid')
env_aniso_norisk_back.showPathData('total-cost',fig,axes,'b','dashed')
env_iso_norisk_go.showPathData('total-cost',fig,axes,'m','solid')
env_iso_norisk_back.showPathData('total-cost',fig,axes,'m','dashed')
env_aniso_wrisk_go.showPathData('total-cost',fig,axes,'c','solid')
env_aniso_wrisk_back.showPathData('total-cost',fig,axes,'c','dashed')
env_iso_wrisk_go.showPathData('total-cost',fig,axes,'r','solid')
env_iso_wrisk_back.showPathData('total-cost',fig,axes,'r','dashed')
plt.style.use('default')


#envBack_bidem.showPathData('total-cost-estimated',fig,axes,'b')
#envGo_biimu.showPathData('total-cost-estimated',fig,axes,'m')
#envBack_biimu.showPathData('total-cost-estimated',fig,axes,'g')

#envGo.show3dDEM()
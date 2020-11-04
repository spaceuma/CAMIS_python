# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:35:43 2020

@author: Richi
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from context import camis
import copy
import zipfile
import yaml
try:
    from scipy import signal
except:
    raise ImportError('ERROR: scipy module could not be imported')

# =============================================================================
## LOADING DEM
# =============================================================================

hiRes_elevationMap = np.loadtxt(\
                        open("data/umaRescueArea/UMARescueArea_1mDEM.csv",\
                             "rb"), delimiter=" ", skiprows=0)

#hiRes_posX = np.loadtxt(open("data/terrainData/UMARescueArea_10cmPosX.csv",\
#                                     "rb"), delimiter=" ", skiprows=0)
#hiRes_posY = np.loadtxt(open("data/terrainData/UMARescueArea_10cmPosY.csv",\
#                                     "rb"), delimiter=" ", skiprows=0)
hiRes = 1.0
offset = np.loadtxt(\
                        open("data/umaRescueArea/UMARescueArea_1mOffset.csv",\
                             "rb"), delimiter=" ", skiprows=0)
print('TEST_DEMO: DEM is loaded')

hexRes = 1.0

env = camis.AnisotropicMap(hiRes_elevationMap, hiRes, hexRes,\
                               offset)

posA = np.asarray([10,180])
posB = np.asarray([180,10])

def computeAllPlannings(anisoMapList):
    anisoMapList[0].executeBiPlanning(posB,posA)
    anisoMapList[1].executeBiPlanning(posA,posB)

def getMapLists(camisInput):
    iso_model = copy.deepcopy(camisInput)
    iso_model.setAsIsotropic()
    env_aniso = [None] * 2
    env_iso = [None] * 2
    env.computeVecCostMap(camisInput)
    enviso = copy.deepcopy(env)
    enviso.computeVecCostMap(iso_model)
    for i in range(2):
        env_aniso[i] = copy.deepcopy(env)
        env_iso[i] = copy.deepcopy(enviso)
    return env_aniso, env_iso



with open("data/sim01/cuadriga_aniso_01.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_01 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD01_scene01, env_isoCUAD01_scene01 = getMapLists(aniso_01)
computeAllPlannings(env_CUAD01_scene01)
computeAllPlannings(env_isoCUAD01_scene01)

with open("data/sim01/cuadriga_aniso_02.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_02 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD02_scene01, env_isoCUAD02_scene01 = getMapLists(aniso_02)
computeAllPlannings(env_CUAD02_scene01)
computeAllPlannings(env_isoCUAD02_scene01)

with open("data/sim01/cuadriga_aniso_03.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_03 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD03_scene01, env_isoCUAD03_scene01 = getMapLists(aniso_03)
computeAllPlannings(env_CUAD03_scene01)
computeAllPlannings(env_isoCUAD03_scene01)

with open("data/sim01/cuadriga_aniso_04.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_04 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD04_scene01, env_isoCUAD04_scene01 = getMapLists(aniso_04)
computeAllPlannings(env_CUAD04_scene01)
computeAllPlannings(env_isoCUAD04_scene01)

with open("data/sim01/cuadriga_aniso_05.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_05 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD05_scene01, env_isoCUAD05_scene01 = getMapLists(aniso_05)
computeAllPlannings(env_CUAD05_scene01)
computeAllPlannings(env_isoCUAD05_scene01)

with open("data/sim01/cuadriga_aniso_06.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_06 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD06_scene01, env_isoCUAD06_scene01 = getMapLists(aniso_06)
computeAllPlannings(env_CUAD06_scene01)
computeAllPlannings(env_isoCUAD06_scene01)

with open("data/sim01/cuadriga_aniso_07.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_07 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD07_scene01, env_isoCUAD07_scene01 = getMapLists(aniso_07)
computeAllPlannings(env_CUAD07_scene01)
computeAllPlannings(env_isoCUAD07_scene01)

with open("data/sim01/cuadriga_aniso_08.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_08 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD08_scene01, env_isoCUAD08_scene01 = getMapLists(aniso_08)
computeAllPlannings(env_CUAD08_scene01)
computeAllPlannings(env_isoCUAD08_scene01)

with open("data/sim01/cuadriga_aniso_09.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_09 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD09_scene01, env_isoCUAD09_scene01 = getMapLists(aniso_09)
computeAllPlannings(env_CUAD09_scene01)
computeAllPlannings(env_isoCUAD09_scene01)


env.show3dDEM()

plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, (ax1,ax2) = plt.subplots(figsize=(5, 4),nrows = 1, ncols = 2, constrained_layout=True)
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*env.hexSlopeMap, 
                 cmap="nipy_spectral",s=4.0)
cbar = fig.colorbar(cc, ax=ax1,shrink=0.6)
cbar.set_label('Steepness (deg)')

def showAnisoPath(mapList, color, ax1, ax2, mode):
    for i,anisomap in enumerate(mapList):
        if i < 1:
            anisomap.showPath(fig,ax1,color,mode)
        else:
            anisomap.showPath(fig,ax2,color,mode)
def showIsoPath(mapList, color, ax1, ax2, mode):
    for i,anisomap in enumerate(mapList):
        if i < 1:
            anisomap.showPath(fig,ax1,color,mode)
        else:
            anisomap.showPath(fig,ax2,color,mode)
            
plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, axes = plt.subplots(figsize=(8, 8), \
      nrows = 1, ncols = 3, \
      sharex = 'all', sharey = 'all')
ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]


plt.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.5, top = 1.0, wspace = 0.01, hspace = 0.05)
fig.text(0.5, 0.005, 'X-axis [m]', ha='center')
fig.text(0.005, 0.5, 'Y-axis [m]', va='center', rotation='vertical')

showAnisoPath(env_CUAD01_scene01, 'r', ax1, ax2,'solid')
showIsoPath(env_isoCUAD01_scene01, 'r', ax3, ax3,'solid')
showAnisoPath(env_CUAD02_scene01, 'c', ax1, ax2,'solid')
showIsoPath(env_isoCUAD02_scene01, 'c', ax3, ax3,'solid')
showAnisoPath(env_CUAD03_scene01, 'lime', ax1, ax2,'solid')
showIsoPath(env_isoCUAD03_scene01, 'lime', ax3, ax3,'solid')
showAnisoPath(env_CUAD04_scene01, 'y', ax1, ax2,'dashed')
showIsoPath(env_isoCUAD04_scene01, 'y', ax3, ax3,'dashed')
showAnisoPath(env_CUAD05_scene01, 'b', ax1, ax2,'dashed')
showIsoPath(env_isoCUAD05_scene01, 'b', ax3, ax3,'dashed')
showAnisoPath(env_CUAD06_scene01, 'g', ax1, ax2,'dashed')
showIsoPath(env_isoCUAD06_scene01, 'g', ax3, ax3,'dashed')
showAnisoPath(env_CUAD07_scene01, 'orange', ax1, ax2,'dotted')
showIsoPath(env_isoCUAD07_scene01, 'orange', ax3, ax3,'dotted')
showAnisoPath(env_CUAD08_scene01, 'm', ax1, ax2,'dotted')
showIsoPath(env_isoCUAD08_scene01, 'm', ax3, ax3,'dotted')   
showAnisoPath(env_CUAD09_scene01, 'k', ax1, ax2,'dotted')
showIsoPath(env_isoCUAD09_scene01, 'k', ax3, ax3,'dotted')
ax1.text(0.5, 0.95, 'Go traverses \n (CAMIS models)', horizontalalignment='center', \
         verticalalignment='center', transform=ax1.transAxes, fontsize = 12,\
         color = 'white')
ax2.text(0.5, 0.95, 'Return traverses \n (CAMIS models)', horizontalalignment='center', \
         verticalalignment='center', transform=ax2.transAxes, fontsize = 12,\
         color = 'white')
ax3.text(0.5, 0.95, 'Go and Return traverses \n (isotropic equivalent models)', horizontalalignment='center', \
         verticalalignment='center', transform=ax3.transAxes, fontsize = 12,\
         color = 'white')
#ax1.legend()
     
for ax in axes:
    cc = ax.scatter(env.hexXmap, env.hexYmap, c = env.hexElevationMap, cmap = cm.gist_earth,s=20)
    ax.scatter(posA[0], posA[1], facecolor = 'r', edgecolor='black', s=60)
    ax.scatter(posB[0], posB[1], facecolor = 'r', edgecolor='black', s=60)
    ax.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
    ax.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
    ax.set_aspect('equal')
fig.tight_layout()



#
#env.executeBiPlanning(posB,posA)
#
#plt.style.use('default')
#plt.rcParams["font.family"] = "Constantia"
#plt.rcParams['mathtext.fontset'] = 'cm'
#plt.rcParams['mathtext.rm'] = 'serif'
#fig, ax = plt.subplots(figsize=(8, 8), \
#      nrows = 1, ncols = 1, \
#      sharex = 'all', sharey = 'all')
#plt.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.5, top = 1.0, wspace = 0.01, hspace = 0.05)
#fig.text(0.5, 0.005, 'X-axis [m]', ha='center')
#fig.text(0.005, 0.5, 'Y-axis [m]', va='center', rotation='vertical')
#
#
#env.showPath(fig,ax,'r','solid')
#
#cc = ax.scatter(env.hexXmap, env.hexYmap, c = env.hexElevationMap, cmap = cm.gist_earth,s=20)
#ax.scatter(posA[0], posA[1], facecolor = 'r', edgecolor='black', s=60)
#ax.scatter(posB[0], posB[1], facecolor = 'r', edgecolor='black', s=60)
#ax.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
#ax.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
#ax.set_aspect('equal')
#fig.tight_layout()















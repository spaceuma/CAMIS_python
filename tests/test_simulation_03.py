# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 11:20:22 2020

@author: rsanchez
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from context import camis
import copy
import zipfile
import yaml

# =============================================================================
## LOADING DEM
# =============================================================================
hiRes_elevationMap = np.loadtxt(\
                        open("data/terrainData/UMATerrainCuesta_10cmDEM.csv",\
                        "rb"), delimiter=" ", skiprows=0)
    
hiRes = 0.1
offset = np.loadtxt(open("data/terrainData/UMATerrainCuesta_10cmOffset.csv",\
                                 "rb"), delimiter=" ", skiprows=0)
env = camis.AnisotropicMap(hiRes_elevationMap[100:400,100:300], hiRes, 0.4,\
                               offset)
env2 = camis.AnisotropicMap(hiRes_elevationMap[400:700,100:300], hiRes, 0.4,\
                               offset)
pos1A = np.asarray([3,27])
pos1B = np.asarray([3,4])
pos1C = np.asarray([17,27])
pos1D = np.asarray([17,4])
print('TEST_DEMO: DEM is loaded')

with open("data/sim01/cuadriga_aniso_01.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_01 = camis.CamisDrivingModel(cuadriga_data)

with open("data/sim01/cuadriga_aniso_02.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_02 = camis.CamisDrivingModel(cuadriga_data)

env.computeVecCostMap(aniso_02)
env2.computeVecCostMap(aniso_02)

env.show3dDEM()
env2.show3dDEM()

plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, (ax1,ax2) = plt.subplots(nrows = 2, ncols = 1, constrained_layout=True)
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*env.hexSlopeMap, cmap="nipy_spectral",s=20)
cc = ax2.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*env2.hexSlopeMap, cmap="nipy_spectral",s=20)
ax1.scatter(pos1A[0], pos1A[1], facecolor = 'r', edgecolor='black', s=60)
ax1.scatter(pos1B[0], pos1B[1], facecolor = 'r', edgecolor='black', s=60)
ax1.scatter(pos1C[0], pos1C[1], facecolor = 'r', edgecolor='black', s=60)
ax1.scatter(pos1D[0], pos1D[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(pos1A[0], pos1A[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(pos1B[0], pos1B[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(pos1C[0], pos1C[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(pos1D[0], pos1D[1], facecolor = 'r', edgecolor='black', s=60)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax1.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
ax1.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
ax2.set_xlim([env2.xMap[0,2], env2.xMap[-1,-4]])
ax2.set_ylim([env2.yMap[0,0], env2.yMap[-1,-1]])
ax1.set_xlabel('X-axis [m]')
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top') 
ax1.set_ylabel('Y-axis [m]')
ax2.set_xlabel('X-axis [m]')
ax2.set_ylabel('Y-axis [m]')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.75, 0.15, 0.05, 0.7])
cbar = fig.colorbar(cc, cax=cbar_ax)
cbar.set_label('Steepness (deg)')
plt.show()
fig, (ax1,ax2) = plt.subplots(nrows = 2, ncols = 1, constrained_layout=True)
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*np.arctan2(env.hexAspectMap[1],env.hexAspectMap[0]), cmap="gist_rainbow",s=20)
cc = ax2.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*np.arctan2(env2.hexAspectMap[1],env2.hexAspectMap[0]), cmap="gist_rainbow",s=20)
ax1.scatter(pos1A[0], pos1A[1], facecolor = 'r', edgecolor='black', s=60)
ax1.scatter(pos1B[0], pos1B[1], facecolor = 'r', edgecolor='black', s=60)
ax1.scatter(pos1C[0], pos1C[1], facecolor = 'r', edgecolor='black', s=60)
ax1.scatter(pos1D[0], pos1D[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(pos1A[0], pos1A[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(pos1B[0], pos1B[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(pos1C[0], pos1C[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(pos1D[0], pos1D[1], facecolor = 'r', edgecolor='black', s=60)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax1.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
ax1.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
ax2.set_xlim([env2.xMap[0,2], env2.xMap[-1,-4]])
ax2.set_ylim([env2.yMap[0,0], env2.yMap[-1,-1]])
ax1.set_xlabel('X-axis [m]')
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top') 
ax1.set_ylabel('Y-axis [m]')
ax2.set_xlabel('X-axis [m]')
ax2.set_ylabel('Y-axis [m]')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.71, 0.15, 0.05, 0.7])
cbar = fig.colorbar(cc, cax=cbar_ax)
cbar.set_label('Aspect (deg)')
plt.show()






def computeAllPlannings(anisoMapList):
    anisoMapList[0].executeBiPlanning(pos1C,pos1A)
    anisoMapList[1].executeBiPlanning(pos1A,pos1C)
    anisoMapList[2].executeBiPlanning(pos1D,pos1A)
    anisoMapList[3].executeBiPlanning(pos1A,pos1D)
    anisoMapList[4].executeBiPlanning(pos1C,pos1B)
    anisoMapList[5].executeBiPlanning(pos1B,pos1C)
    anisoMapList[6].executeBiPlanning(pos1D,pos1B)
    anisoMapList[7].executeBiPlanning(pos1B,pos1D) 

def getMapLists(camisInput):
    iso_model = copy.deepcopy(camisInput)
    iso_model.setAsIsotropic()
    env_aniso = [None] * 8
    env_iso = [None] * 8
    env.computeVecCostMap(camisInput)
    enviso = copy.deepcopy(env)
    enviso.computeVecCostMap(iso_model)
    for i in range(8):
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

plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, (axes1,axes2,axes3) = plt.subplots(nrows = 3, ncols = 2)
ax1 = axes1[0]
ax2 = axes1[1]
ax3 = axes2[0]
ax4 = axes2[1]
ax5 = axes3[0]
ax6 = axes3[1]
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = env.hexElevationMap, cmap=cm.gist_earth,s=20)
cc = ax2.scatter(env.hexXmap, env.hexYmap, c = env2.hexElevationMap, cmap=cm.gist_earth,s=20)
cc = ax3.scatter(env.hexXmap, env.hexYmap, c = env.hexElevationMap, cmap=cm.gist_earth,s=20)
cc = ax4.scatter(env.hexXmap, env.hexYmap, c = env2.hexElevationMap, cmap=cm.gist_earth,s=20)
cc = ax5.scatter(env.hexXmap, env.hexYmap, c = env.hexElevationMap, cmap=cm.gist_earth,s=20)
cc = ax6.scatter(env.hexXmap, env.hexYmap, c = env2.hexElevationMap, cmap=cm.gist_earth,s=20)
for i,anisomap in enumerate(env_CUAD01_scene01):
    if i%2==0:
        if i < 4:
            anisomap.showPath(fig,ax1,'r','solid')
        else:
            anisomap.showPath(fig,ax3,'r','solid')
    else:
        if i < 4:
            anisomap.showPath(fig,ax1,'r','dashed')
        else:
            anisomap.showPath(fig,ax3,'r','dashed')
for i,anisomap in enumerate(env_CUAD02_scene01):
    if i%2==0:
        if i < 4:
            anisomap.showPath(fig,ax1,'c','solid')
        else:
            anisomap.showPath(fig,ax3,'c','solid')
    else:
        if i < 4:
            anisomap.showPath(fig,ax1,'c','dashed')
        else:
            anisomap.showPath(fig,ax3,'c','dashed')
for i,anisomap in enumerate(env_CUAD03_scene01):
    if i%2==0:
        if i < 4:
            anisomap.showPath(fig,ax1,'g','solid')
        else:
            anisomap.showPath(fig,ax3,'g','solid')
    else:
        if i < 4:
            anisomap.showPath(fig,ax1,'g','dashed')
        else:
            anisomap.showPath(fig,ax3,'g','dashed')
for isomap in env_isoCUAD01_scene01:
    isomap.showPath(fig,ax5,'r','solid')
for isomap in env_isoCUAD02_scene01:
    isomap.showPath(fig,ax5,'c','solid')
for isomap in env_isoCUAD03_scene01:
    isomap.showPath(fig,ax5,'g','solid')
ax1.scatter(pos1A[0], pos1A[1], facecolor = 'r', edgecolor='black', s=60)
#ax1.scatter(pos1B[0], pos1B[1], facecolor = 'r', edgecolor='black', s=60)
ax1.scatter(pos1C[0], pos1C[1], facecolor = 'r', edgecolor='black', s=60)
ax1.scatter(pos1D[0], pos1D[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(pos1A[0], pos1A[1], facecolor = 'r', edgecolor='black', s=60)
#ax2.scatter(pos1B[0], pos1B[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(pos1C[0], pos1C[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(pos1D[0], pos1D[1], facecolor = 'r', edgecolor='black', s=60)
#ax3.scatter(pos1A[0], pos1A[1], facecolor = 'r', edgecolor='black', s=60)
ax3.scatter(pos1B[0], pos1B[1], facecolor = 'r', edgecolor='black', s=60)
ax3.scatter(pos1C[0], pos1C[1], facecolor = 'r', edgecolor='black', s=60)
ax3.scatter(pos1D[0], pos1D[1], facecolor = 'r', edgecolor='black', s=60)
#ax4.scatter(pos1A[0], pos1A[1], facecolor = 'r', edgecolor='black', s=60)
ax4.scatter(pos1B[0], pos1B[1], facecolor = 'r', edgecolor='black', s=60)
ax4.scatter(pos1C[0], pos1C[1], facecolor = 'r', edgecolor='black', s=60)
ax4.scatter(pos1D[0], pos1D[1], facecolor = 'r', edgecolor='black', s=60)
ax5.scatter(pos1A[0], pos1A[1], facecolor = 'r', edgecolor='black', s=60)
ax5.scatter(pos1B[0], pos1B[1], facecolor = 'r', edgecolor='black', s=60)
ax5.scatter(pos1C[0], pos1C[1], facecolor = 'r', edgecolor='black', s=60)
ax5.scatter(pos1D[0], pos1D[1], facecolor = 'r', edgecolor='black', s=60)
ax6.scatter(pos1A[0], pos1A[1], facecolor = 'r', edgecolor='black', s=60)
ax6.scatter(pos1B[0], pos1B[1], facecolor = 'r', edgecolor='black', s=60)
ax6.scatter(pos1C[0], pos1C[1], facecolor = 'r', edgecolor='black', s=60)
ax6.scatter(pos1D[0], pos1D[1], facecolor = 'r', edgecolor='black', s=60)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax4.set_aspect('equal')
ax5.set_aspect('equal')
ax6.set_aspect('equal')
ax1.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
ax1.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
ax2.set_xlim([env2.xMap[0,2], env2.xMap[-1,-4]])
ax2.set_ylim([env2.yMap[0,0], env2.yMap[-1,-1]])
ax3.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
ax3.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
ax4.set_xlim([env2.xMap[0,2], env2.xMap[-1,-4]])
ax4.set_ylim([env2.yMap[0,0], env2.yMap[-1,-1]])
ax5.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
ax5.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
ax6.set_xlim([env2.xMap[0,2], env2.xMap[-1,-4]])
ax6.set_ylim([env2.yMap[0,0], env2.yMap[-1,-1]])
ax1.set_xlabel('X-axis [m]')
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top') 
ax1.set_ylabel('Y-axis [m]')
ax2.set_xlabel('X-axis [m]')
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position('top') 
ax3.set_ylabel('Y-axis [m]')
ax3.set_xlabel('X-axis [m]')
ax5.set_xlabel('X-axis [m]')
ax5.set_ylabel('Y-axis [m]')
ax6.set_xlabel('X-axis [m]')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.8, 0.15, 0.05, 0.7])
cbar = fig.colorbar(cc, cax=cbar_ax)
cbar.set_label('Elevation [m]')
plt.show()
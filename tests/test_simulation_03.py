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
from scipy.optimize import curve_fit

# =============================================================================
## LOADING DEM
# =============================================================================
hiRes_elevationMap = np.loadtxt(\
                        open("data/terrainData/UMATerrainCuesta_10cmDEM.csv",\
                        "rb"), delimiter=" ", skiprows=0)
    
hiRes = 0.1
offset = np.loadtxt(open("data/terrainData/UMATerrainCuesta_10cmOffset.csv",\
                                 "rb"), delimiter=" ", skiprows=0)
env = camis.AnisotropicMap(hiRes_elevationMap[50:850,80:325], hiRes, 0.4,\
                               offset)
posA = np.asarray([5,70])
posB = np.asarray([20,55])
posC = np.asarray([5,40])
posD = np.asarray([20,25])
posE = np.asarray([5,10])
print('TEST_DEMO: DEM is loaded')

with open("data/sim01/cuadriga_aniso_01.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_01 = camis.CamisDrivingModel(cuadriga_data)

with open("data/sim01/cuadriga_aniso_02.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_02 = camis.CamisDrivingModel(cuadriga_data)

env.computeVecCostMap(aniso_02)

env.show3dDEM()

def showPoints(ax):
    ax.annotate('A',
                xy=(posA[0], posA[1]),
                xytext=(-4, 0),  # 3 points vertical offset
                textcoords="offset points",
                ha='right', va='bottom',color='white')
    ax.annotate('B',
                    xy=(posB[0], posB[1]),
                    xytext=(4, 0),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='left', va='bottom',color='white')
    ax.annotate('C',
                    xy=(posC[0], posC[1]),
                    xytext=(-4, 0),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='right', va='bottom',color='white')
    ax.annotate('D',
                    xy=(posD[0], posD[1]),
                    xytext=(4, 0),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='left', va='bottom',color='white')
    ax.annotate('E',
                    xy=(posE[0], posE[1]),
                    xytext=(-4, 0),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='right', va='bottom',color='white')

plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, (ax1,ax2) = plt.subplots(figsize=(5, 4),nrows = 1, ncols = 2, constrained_layout=True)
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*env.hexSlopeMap, 
                 cmap="nipy_spectral",s=.2)
ax1.scatter(posA[0], posA[1], facecolor = 'r', edgecolor='black', s=60)
ax1.scatter(posB[0], posB[1], facecolor = 'r', edgecolor='black', s=60)
ax1.scatter(posC[0], posC[1], facecolor = 'r', edgecolor='black', s=60)
ax1.scatter(posD[0], posD[1], facecolor = 'r', edgecolor='black', s=60)
ax1.scatter(posE[0], posE[1], facecolor = 'r', edgecolor='black', s=60)
showPoints(ax1)
ax1.set_aspect('equal')
ax1.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
ax1.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
ax1.set_xlabel('X-axis [m]')
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top') 
ax1.set_ylabel('Y-axis [m]')
cbar = fig.colorbar(cc, ax=ax1,shrink=0.6)
cbar.set_label('Steepness (deg)')
cc = ax2.scatter(env.hexXmap, env.hexYmap, 
                 c = 180/np.pi*np.arctan2(env.hexAspectMap[1],env.hexAspectMap[0]),
                 cmap="gist_rainbow",s=.2)
ax2.scatter(posA[0], posA[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(posB[0], posB[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(posC[0], posC[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(posD[0], posD[1], facecolor = 'r', edgecolor='black', s=60)
ax2.scatter(posE[0], posE[1], facecolor = 'r', edgecolor='black', s=60)
showPoints(ax2)
ax2.set_aspect('equal')
ax2.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
ax2.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
ax2.set_xlabel('X-axis [m]')
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position('top') 
ax2.set_ylabel('Y-axis [m]')
cbar = fig.colorbar(cc, ax=ax2,shrink=0.6)
cbar.set_label('Aspect (deg)')
plt.show()


################# SLIP FUNCTIONS

steepnessArray = np.arange(0, 31.5, 0.5)
s1 = 4*np.sin(steepnessArray/31.0*45*np.pi/180.0)**4
s2 = 1.95*np.sin(steepnessArray/31.0*45*np.pi/180.0)**2
s3 = 0.95*np.sin(steepnessArray/31.0*90*np.pi/180.0)
fig, ax = plt.subplots()
ax.plot(steepnessArray, s1, label = 's1')
ax.plot(steepnessArray, s2, label = 's2')
ax.plot(steepnessArray, s3, label = 's3')
ax.legend()

def func(x,a,b,c):
    return a*x**3 + b*x**2 + c*x

roots1,_ = curve_fit(func,steepnessArray, s1)
roots2,_  = curve_fit(func,steepnessArray, s2)
roots3,_  = curve_fit(func,steepnessArray, s3)

roots1 = np.insert(roots1,3,0.0)
roots2 = np.insert(roots2,3,0.0)
roots3 = np.insert(roots3,3,0.0)

fig, ax = plt.subplots()
ax.plot(steepnessArray, np.polyval(roots1,steepnessArray), label = 's1')
ax.plot(steepnessArray, np.polyval(roots2,steepnessArray), label = 's2')
ax.plot(steepnessArray, np.polyval(roots3,steepnessArray), label = 's3')
ax.legend()


######## MODELS ####

def computeAllPlannings(anisoMapList):
    anisoMapList[0].executeBiPlanning(posB,posA)
    anisoMapList[1].executeBiPlanning(posC,posB)
    anisoMapList[2].executeBiPlanning(posD,posC)
    anisoMapList[3].executeBiPlanning(posE,posD)
    anisoMapList[4].executeBiPlanning(posD,posE)
    anisoMapList[5].executeBiPlanning(posC,posD)
    anisoMapList[6].executeBiPlanning(posB,posC)
    anisoMapList[7].executeBiPlanning(posA,posB)

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

with open("data/sim01/cuadriga_aniso_04.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_04 = camis.CamisDrivingModel(cuadriga_data)
aniso_04.showDirCosts()
aniso_04.showAnisotropy()
env_CUAD04_scene01, env_isoCUAD04_scene01 = getMapLists(aniso_04)
computeAllPlannings(env_CUAD04_scene01)
computeAllPlannings(env_isoCUAD04_scene01)

with open("data/sim01/cuadriga_aniso_05.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_05 = camis.CamisDrivingModel(cuadriga_data)
aniso_05.showDirCosts()
aniso_05.showAnisotropy()
env_CUAD05_scene01, env_isoCUAD05_scene01 = getMapLists(aniso_05)
computeAllPlannings(env_CUAD05_scene01)
computeAllPlannings(env_isoCUAD05_scene01)

with open("data/sim01/cuadriga_aniso_06.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_06 = camis.CamisDrivingModel(cuadriga_data)
aniso_06.showDirCosts()
aniso_06.showAnisotropy()
env_CUAD06_scene01, env_isoCUAD06_scene01 = getMapLists(aniso_06)
computeAllPlannings(env_CUAD06_scene01)
computeAllPlannings(env_isoCUAD06_scene01)

with open("data/sim01/cuadriga_aniso_07.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_07 = camis.CamisDrivingModel(cuadriga_data)
aniso_07.showDirCosts()
aniso_07.showAnisotropy()
env_CUAD07_scene01, env_isoCUAD07_scene01 = getMapLists(aniso_07)
computeAllPlannings(env_CUAD07_scene01)
computeAllPlannings(env_isoCUAD07_scene01)

with open("data/sim01/cuadriga_aniso_08.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_08 = camis.CamisDrivingModel(cuadriga_data)
aniso_08.showDirCosts()
aniso_08.showAnisotropy()
env_CUAD08_scene01, env_isoCUAD08_scene01 = getMapLists(aniso_08)
computeAllPlannings(env_CUAD08_scene01)
computeAllPlannings(env_isoCUAD08_scene01)

with open("data/sim01/cuadriga_aniso_09.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_09 = camis.CamisDrivingModel(cuadriga_data)
aniso_09.showDirCosts()
aniso_09.showAnisotropy()
env_CUAD09_scene01, env_isoCUAD09_scene01 = getMapLists(aniso_09)
computeAllPlannings(env_CUAD09_scene01)
computeAllPlannings(env_isoCUAD09_scene01)


def showAnisoPath(mapList, color, ax1, ax2, mode):
    for i,anisomap in enumerate(mapList):
        if i < 4:
            anisomap.showPath(fig,ax1,color,mode)
        else:
            anisomap.showPath(fig,ax2,color,mode)
def showIsoPath(mapList, color, ax1, ax2, mode):
    for i,anisomap in enumerate(mapList):
        if i < 4:
            anisomap.showPath(fig,ax1,color,mode)
        else:
            anisomap.showPath(fig,ax2,color,mode)

######################## SHOW PATHS ##################
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

showAnisoPath(env_CUAD01_scene01, 'r', ax1, ax2)
showIsoPath(env_isoCUAD01_scene01, 'r', ax3, ax3)
showAnisoPath(env_CUAD02_scene01, 'c', ax1, ax2)
showIsoPath(env_isoCUAD02_scene01, 'c', ax3, ax3)
showAnisoPath(env_CUAD03_scene01, 'lime', ax1, ax2)
showIsoPath(env_isoCUAD03_scene01, 'lime', ax3, ax3)
showAnisoPath(env_CUAD04_scene01, 'y', ax1, ax2)
showIsoPath(env_isoCUAD04_scene01, 'y', ax3, ax3)
showAnisoPath(env_CUAD05_scene01, 'b', ax1, ax2)
showIsoPath(env_isoCUAD05_scene01, 'b', ax3, ax3)
showAnisoPath(env_CUAD06_scene01, 'g', ax1, ax2)
showIsoPath(env_isoCUAD06_scene01, 'g', ax3, ax3)
showAnisoPath(env_CUAD07_scene01, 'orange', ax1, ax2)
showIsoPath(env_isoCUAD07_scene01, 'orange', ax3, ax3)
showAnisoPath(env_CUAD08_scene01, 'm', ax1, ax2)
showIsoPath(env_isoCUAD08_scene01, 'm', ax3, ax3)   
showAnisoPath(env_CUAD09_scene01, 'k', ax1, ax2)
showIsoPath(env_isoCUAD09_scene01, 'k', ax3, ax3)  
showPoints(ax1)
showPoints(ax2)
showPoints(ax3)
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
    ax.scatter(posC[0], posC[1], facecolor = 'r', edgecolor='black', s=60)
    ax.scatter(posD[0], posD[1], facecolor = 'r', edgecolor='black', s=60)
    ax.scatter(posE[0], posE[1], facecolor = 'r', edgecolor='black', s=60)
    ax.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
    ax.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
    ax.set_aspect('equal')
fig.tight_layout()

#fig.subplots_adjust(right=0.87)
cbar_ax = fig.add_axes([0.1, 0.12, 0.2, 0.02])
legend_ax = fig.add_axes([0.5, 0.11, 0.4, 0.05])
legend_ax.plot(0,0,'r',label='CUAD01')
legend_ax.plot(0,0,'c',label='CORI01')
legend_ax.plot(0,0,'lime',label='CORI02')
legend_ax.plot(0,0,'y',label='CSLI01')
legend_ax.plot(0,0,'b',label='CSLI02')
legend_ax.plot(0,0,'g',label='CSLI03')
legend_ax.plot(0,0,'orange',label='CTIM01')
legend_ax.plot(0,0,'m',label='CTIM02')
legend_ax.plot(0,0,'k',label='CTIM03')
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(b=False)
plt.axis('off')
cbar = fig.colorbar(cc, cax=cbar_ax, orientation = 'horizontal')
cbar.set_label('Elevation [m]')
legend_ax.legend(ncol=3)
plt.show()



############### ORIENTATION ###############

def showOrientation(env,ax,color,name):
    for i in range(8):
        if i == 0:
            offSet = 0
        else:
            offSet = offSet + env[i-1].pathTravDist[-1]
        ax.fill_between(offSet + np.asarray(env[i].pathTravDist), 0, \
                  180.0/np.pi*np.asarray(env[i].pathSlope), facecolor = 'y',\
                  alpha = 0.5)
        ax.plot(offSet + np.asarray(env[i].pathTravDist), \
                      180.0/np.pi*np.asarray(env[i].pathSlope), 'y')
        ax.fill_between(offSet + np.asarray(env[i].pathTravDist), 0, \
                  -180.0/np.pi*np.asarray(env[i].pathSlope), facecolor = 'y',\
                  alpha = 0.5)
        ax.plot(offSet + np.asarray(env[i].pathTravDist), \
                      -180.0/np.pi*np.asarray(env[i].pathSlope), 'y')
        ax.fill_between(offSet + np.asarray(env[i].pathTravDist), 0, \
                  180.0/np.pi*np.asarray(env[i].pathPitch), facecolor = 'c',\
                  alpha = 0.5)
        ax.plot(offSet + np.asarray(env[i].pathTravDist), \
                      180.0/np.pi*np.asarray(env[i].pathPitch), 'c')
        ax.fill_between(offSet + np.asarray(env[i].pathTravDist), 0, \
                  180.0/np.pi*env[i].pathRoll*np.sign(env[i].pathBeta), facecolor = 'r',\
                  alpha = 0.5)
        ax.plot(offSet + np.asarray(env[i].pathTravDist), \
                      180.0/np.pi*env[i].pathRoll*np.sign(env[i].pathBeta), 'r')
    ax.annotate(name,
                xy=(230, 5),
                xytext=(-4, 0),  # 3 points vertical offset
                textcoords="offset points",
                ha='right', va='bottom',color=color)

plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, rowaxes = plt.subplots(figsize=(7, 7), nrows = 6, ncols = 1, \
     sharex = 'all', sharey = 'all')

rowaxes[5].set_xlabel('Traversed Distance [m]')

plt.subplots_adjust(left = 0.085, right = 0.95, bottom = 0.075, top = 0.9, wspace = 0.1, hspace = 0.075)
fig.text(0.02, 0.5, 'Orientation Angle [deg]', va='center', rotation='vertical')
rowaxes[0].set_xlim([0,230])

showOrientation(env_CUAD01_scene01, rowaxes[0], 'r', 'CUAD01')
showOrientation(env_CUAD02_scene01, rowaxes[2], 'c', 'CORI01')
showOrientation(env_CUAD04_scene01, rowaxes[4], 'y', 'CORI02')

showOrientation(env_isoCUAD01_scene01, rowaxes[1], 'r', 'CUAD01 (isotropic)')
showOrientation(env_isoCUAD02_scene01, rowaxes[3], 'c', 'CORI01 (isotropic)')
showOrientation(env_isoCUAD04_scene01, rowaxes[5], 'y', 'CORI02 (isotropic)')

plt.grid(b=True, which = 'minor')
legend_ax = fig.add_axes([0.2, 0.92, 0.6, 0.07])
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(b=False)
xArray = np.arange(0, 2*np.pi, 0.01)
legend_ax.plot(xArray, np.sin(xArray),'y')
legend_ax.plot(xArray, -np.sin(xArray),'y')
legend_ax.fill_between(xArray, -np.sin(xArray), \
                  np.sin(xArray), facecolor = 'y',\
                  alpha = 0.5)
legend_ax.plot(xArray+2*np.pi, np.sin(xArray),'c')
legend_ax.fill_between(xArray+2*np.pi, 0, \
                  np.sin(xArray), facecolor = 'c',\
                  alpha = 0.5)
legend_ax.plot(xArray+4*np.pi, np.sin(xArray),'r')
legend_ax.fill_between(xArray+4*np.pi, 0, \
                  np.sin(xArray), facecolor = 'r',\
                  alpha = 0.5)
legend_ax.text(np.pi/2, 0.0, 'Steepness', horizontalalignment='center', \
         verticalalignment='center', fontsize = 10,\
         color = 'white')
legend_ax.text(3*np.pi/2, -0.5, '-α', horizontalalignment='center', \
         verticalalignment='center', fontsize = 12,\
         color = 'white')
legend_ax.text(3*np.pi/2, 0.5, 'α', horizontalalignment='center', \
         verticalalignment='center', fontsize = 12,\
         color = 'white')
legend_ax.text(5*np.pi/2, 0.5, 'Pitch', horizontalalignment='center', \
         verticalalignment='center', fontsize = 10,\
         color = 'white')
legend_ax.text(7*np.pi/2, -0.5, 'θ', horizontalalignment='center', \
         verticalalignment='center', fontsize = 10,\
         color = 'white')
legend_ax.text(9*np.pi/2, 0.5, 'Roll', horizontalalignment='center', \
         verticalalignment='center', fontsize = 10,\
         color = 'white')
legend_ax.text(11*np.pi/2, -0.5, 'Φ', horizontalalignment='center', \
         verticalalignment='center', fontsize = 10,\
         color = 'white')
legend_ax.set_facecolor('w')
plt.show()




################ ENERGY PLOT ##################

coeffsLabels = ['CUAD01','CUAD02','CUAD03','CSLI01','CSLI02','CSLI03']
anisoTotalCost = [0,0,0,0,0,0]
isoTotalCost = [0,0,0,0,0,0]
for i in range(8):
    anisoTotalCost[0] = anisoTotalCost[0] + env_CUAD01_scene01[i].pathComputedTotalCost[-1]
    anisoTotalCost[1] = anisoTotalCost[1] + env_CUAD02_scene01[i].pathComputedTotalCost[-1]
    anisoTotalCost[2] = anisoTotalCost[2] + env_CUAD03_scene01[i].pathComputedTotalCost[-1]
    anisoTotalCost[3] = anisoTotalCost[3] + env_CUAD04_scene01[i].pathComputedTotalCost[-1]
    anisoTotalCost[4] = anisoTotalCost[4] + env_CUAD05_scene01[i].pathComputedTotalCost[-1]
    anisoTotalCost[5] = anisoTotalCost[5] + env_CUAD06_scene01[i].pathComputedTotalCost[-1]
    isoTotalCost[0] = isoTotalCost[0] + env_isoCUAD01_scene01[i].pathComputedTotalCost[-1]
    isoTotalCost[1] = isoTotalCost[1] + env_isoCUAD02_scene01[i].pathComputedTotalCost[-1]
    isoTotalCost[2] = isoTotalCost[2] + env_isoCUAD03_scene01[i].pathComputedTotalCost[-1]
    isoTotalCost[3] = isoTotalCost[3] + env_isoCUAD04_scene01[i].pathComputedTotalCost[-1]
    isoTotalCost[4] = isoTotalCost[4] + env_isoCUAD05_scene01[i].pathComputedTotalCost[-1]
    isoTotalCost[5] = isoTotalCost[5] + env_isoCUAD06_scene01[i].pathComputedTotalCost[-1]


x = np.arange(len(coeffsLabels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
ax.bar(x - width/2, anisoTotalCost, width, label='Isotropic (ρ = 0.8)', color='r')
ax.bar(x + width/2, isoTotalCost, width, label='Isotropic (ρ = 0.8)', color = 'b')
#autolabel(rects1)
#autolabel(rects2)
ax.set_ylabel('Total Cost [As]')
ax.set_xlabel('CAMIS Model')
ax.set_xticks(x)
ax.set_xticklabels(coeffsLabels)
#ax.legend()
plt.show()

################ ELAPSED TIME ##################

coeffsLabels = ['CUAD01','CUAD02','CUAD03','CSLI01','CSLI02','CSLI03','CTIM01','CTIM02','CTIM03']
anisoTime = [0,0,0,0,0,0,0,0,0]
isoTime = [0,0,0,0,0,0,0,0,0]
for i in range(8):
    anisoTime[0] = anisoTime[0] + env_CUAD01_scene01[i].elapsedTime
    anisoTime[1] = anisoTime[1] + env_CUAD02_scene01[i].elapsedTime
    anisoTime[2] = anisoTime[2] + env_CUAD03_scene01[i].elapsedTime
    anisoTime[3] = anisoTime[3] + env_CUAD04_scene01[i].elapsedTime
    anisoTime[4] = anisoTime[4] + env_CUAD05_scene01[i].elapsedTime
    anisoTime[5] = anisoTime[5] + env_CUAD06_scene01[i].elapsedTime
    anisoTime[6] = anisoTime[6] + env_CUAD07_scene01[i].elapsedTime
    anisoTime[7] = anisoTime[7] + env_CUAD08_scene01[i].elapsedTime
    anisoTime[8] = anisoTime[8] + env_CUAD09_scene01[i].elapsedTime
    isoTime[0] = isoTime[0] + env_isoCUAD01_scene01[i].elapsedTime
    isoTime[1] = isoTime[1] + env_isoCUAD02_scene01[i].elapsedTime
    isoTime[2] = isoTime[2] + env_isoCUAD03_scene01[i].elapsedTime
    isoTime[3] = isoTime[3] + env_isoCUAD04_scene01[i].elapsedTime
    isoTime[4] = isoTime[4] + env_isoCUAD05_scene01[i].elapsedTime
    isoTime[5] = isoTime[5] + env_isoCUAD06_scene01[i].elapsedTime
    isoTime[6] = isoTime[6] + env_isoCUAD07_scene01[i].elapsedTime
    isoTime[7] = isoTime[7] + env_isoCUAD08_scene01[i].elapsedTime
    isoTime[8] = isoTime[8] + env_isoCUAD09_scene01[i].elapsedTime


x = np.arange(len(coeffsLabels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
ax.bar(x - width/2, anisoTime, width, label='Isotropic (ρ = 0.8)', color='r')
ax.bar(x + width/2, isoTime, width, label='Isotropic (ρ = 0.8)', color = 'b')
#autolabel(rects1)
#autolabel(rects2)
ax.set_ylabel('Total Cost [As]')
ax.set_xlabel('CAMIS Model')
ax.set_xticks(x)
ax.set_xticklabels(coeffsLabels)
#ax.legend()
plt.show()

################ TIME ##################

coeffsLabels = ['CTIM01','CTIM02','CTIM03']
anisoTime = [0,0,0]
isoTime = [0,0,0]
for i in range(8):
    anisoTime[0] = anisoTime[0] + env_CUAD07_scene01[i].pathComputedTotalCost[-1]
    anisoTime[1] = anisoTime[1] + env_CUAD08_scene01[i].pathComputedTotalCost[-1]
    anisoTime[2] = anisoTime[2] + env_CUAD09_scene01[i].pathComputedTotalCost[-1]
    isoTime[0] = isoTime[0] + env_isoCUAD07_scene01[i].pathComputedTotalCost[-1]
    isoTime[1] = isoTime[1] + env_isoCUAD08_scene01[i].pathComputedTotalCost[-1]
    isoTime[2] = isoTime[2] + env_isoCUAD09_scene01[i].pathComputedTotalCost[-1]


x = np.arange(len(coeffsLabels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
ax.bar(x - width/2, anisoTime, width, label='Isotropic (ρ = 0.8)', color='r')
ax.bar(x + width/2, isoTime, width, label='Isotropic (ρ = 0.8)', color = 'b')
#autolabel(rects1)
#autolabel(rects2)
ax.set_ylabel('Total Cost [As]')
ax.set_xlabel('CAMIS Model')
ax.set_xticks(x)
ax.set_xticklabels(coeffsLabels)
#ax.legend()
plt.show()







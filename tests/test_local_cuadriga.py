# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 12:32:04 2020

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
from mpl_toolkits.mplot3d import Axes3D

posA = np.asarray([7,23])
posC = np.asarray([7,5])
posB = np.asarray([22,22])
posD = np.asarray([22,7])


#posD = np.asarray([20,10])

# =============================================================================
## LOADING DEM
# =============================================================================
hiRes_elevationMap = np.loadtxt(\
                        open("data/terrainData/UMATerrainCuesta_10cmDEM.csv",\
                        "rb"), delimiter=" ", skiprows=0)
    
hiRes = 0.1
offset = np.loadtxt(open("data/terrainData/UMATerrainCuesta_10cmOffset.csv",\
                                 "rb"), delimiter=" ", skiprows=0)



hexRes = 1.0/np.sqrt(6)

localDEM = hiRes_elevationMap[260:510,350:600]
offset = offset + [350*0.1, 260*0.1]
occupancy_radius = 0.5
tracking_error = 0.5
env = camis.AnisotropicMap(hiRes_elevationMap[260:510,350:600], hiRes, hexRes,\
                               offset, occupancy_radius, tracking_error)

#offset = offset + [350*0.1, 650*0.1]
#env = camis.AnisotropicMap(hiRes_elevationMap[650:1000,350:600], hiRes, hexRes,\
#                               offset)

def computeAllPlannings(anisoMapList):
    anisoMapList[0].executeBiPlanning(posB,posA)
    anisoMapList[1].executeBiPlanning(posC,posB)
    anisoMapList[2].executeBiPlanning(posD,posC)
    anisoMapList[3].executeBiPlanning(posA,posD)
    anisoMapList[4].executeBiPlanning(posD,posA)
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

with open("data/localCuadrigaTest/cuadriga_aniso_01.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_01 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD01_scene01, env_isoCUAD01_scene01 = getMapLists(aniso_01)
fig, ax = plt.subplots(constrained_layout=True)
#env_CUAD01_scene01[0].showMap('conv-slope',fig,ax)
cc = ax.scatter(env.hexXmap, env.hexYmap, c = env_CUAD01_scene01[0].hexSlopeMap, cmap = cm.gist_earth,s=60)
ax.scatter(posA[0], posA[1], facecolor = 'r', edgecolor='black', s=60)
ax.text(posA[0]-1, posA[1],'A',color = 'white')
ax.scatter(posB[0], posB[1], facecolor = 'r', edgecolor='black', s=60)
ax.text(posB[0]+1, posB[1],'B',color = 'white')
ax.scatter(posC[0], posC[1], facecolor = 'r', edgecolor='black', s=60)
ax.text(posC[0]-1, posC[1],'C',color = 'white')
ax.scatter(posD[0], posD[1], facecolor = 'r', edgecolor='black', s=60)
ax.text(posD[0]+1, posD[1],'D',color = 'white')
ax.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
ax.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
ax.set_aspect('equal')
computeAllPlannings(env_CUAD01_scene01)
computeAllPlannings(env_isoCUAD01_scene01)

print('TEST: Finished CUAD01 plannings')

with open("data/localCuadrigaTest/cuadriga_aniso_02.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_02 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD02_scene01, env_isoCUAD02_scene01 = getMapLists(aniso_02)
computeAllPlannings(env_CUAD02_scene01)
computeAllPlannings(env_isoCUAD02_scene01)

print('TEST: Finished CUAD02 plannings')

with open("data/localCuadrigaTest/cuadriga_aniso_03.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_03 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD03_scene01, env_isoCUAD03_scene01 = getMapLists(aniso_03)
computeAllPlannings(env_CUAD03_scene01)
computeAllPlannings(env_isoCUAD03_scene01)

print('TEST: Finished CUAD03 plannings')

with open("data/localCuadrigaTest/cuadriga_aniso_04.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_04 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD04_scene01, env_isoCUAD04_scene01 = getMapLists(aniso_04)
computeAllPlannings(env_CUAD04_scene01)
computeAllPlannings(env_isoCUAD04_scene01)

print('TEST: Finished CUAD04 plannings')

with open("data/localCuadrigaTest/cuadriga_aniso_05.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_05 = camis.CamisDrivingModel(cuadriga_data)
env_CUAD05_scene01, env_isoCUAD05_scene01 = getMapLists(aniso_05)
computeAllPlannings(env_CUAD05_scene01)
computeAllPlannings(env_isoCUAD05_scene01)

print('TEST: Finished CUAD05 plannings')

#with open("data/localCuadrigaTest/cuadriga_aniso_06.yml", 'r') as file:
#    cuadriga_data = yaml.full_load(file)
#aniso_06 = camis.CamisDrivingModel(cuadriga_data)
#env_CUAD06_scene01, env_isoCUAD06_scene01 = getMapLists(aniso_06)
#computeAllPlannings(env_CUAD06_scene01)
#computeAllPlannings(env_isoCUAD06_scene01)
#
#print('TEST: Finished CUAD06 plannings')
#
#
#with open("data/localCuadrigaTest/cuadriga_aniso_07.yml", 'r') as file:
#    cuadriga_data = yaml.full_load(file)
#aniso_07 = camis.CamisDrivingModel(cuadriga_data)
#env_CUAD07_scene01, env_isoCUAD07_scene01 = getMapLists(aniso_07)
#computeAllPlannings(env_CUAD07_scene01)
#computeAllPlannings(env_isoCUAD07_scene01)
#
#print('TEST: Finished CUAD07 plannings')

#env.show3dDEM()

plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, ax1 = plt.subplots(figsize=(5, 4),nrows = 1, ncols = 1, constrained_layout=True)
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*env.hexSlopeMap, 
                 cmap="nipy_spectral",s=16.0)
cbar = fig.colorbar(cc, ax=ax1,shrink=0.6)
cc.set_clim(0,30.0)
cbar.set_label('Steepness (deg)')





######################## SHOW PATHS ##################
def showAnisoPath(mapList, color, ax1, ax2, mode):
    for i,anisomap in enumerate(mapList):
        if i < 4:
            anisomap.showPath(fig,ax1,color,mode)
        else:
            anisomap.showPath(fig,ax2,color,mode)
def showIsoPath(mapList, color, ax1, mode):
    for i,isomap in enumerate(mapList):
        if i < 4:
            isomap.showPath(fig,ax1,color,mode)
            
plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, axes = plt.subplots(figsize=(12, 4.5), \
      nrows = 1, ncols = 3, \
      sharex = 'all', sharey = 'all')
ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]

fig.text(0.5, 0.005, 'X-axis [m]', ha='center')
fig.text(0.005, 0.5, 'Y-axis [m]', va='center', rotation='vertical')

showAnisoPath(env_CUAD01_scene01, 'r', ax1, ax2, 'solid')
showAnisoPath(env_CUAD02_scene01, 'm', ax1, ax2, 'solid')
showAnisoPath(env_CUAD03_scene01, 'b', ax1, ax2, 'solid')
showAnisoPath(env_CUAD04_scene01, 'c', ax1, ax2, 'solid')
showAnisoPath(env_CUAD05_scene01, 'g', ax1, ax2, 'solid')
showIsoPath(env_isoCUAD01_scene01, 'r', ax3, 'solid')
showIsoPath(env_isoCUAD02_scene01, 'm', ax3, 'solid')
showIsoPath(env_isoCUAD03_scene01, 'b', ax3, 'solid')
showIsoPath(env_isoCUAD04_scene01, 'c', ax3, 'solid')
showIsoPath(env_isoCUAD05_scene01, 'g', ax3, 'solid')
#ax1.text(0.01, 0.1, 'Anisotropic A-B-C-D-A', horizontalalignment='left', \
#         verticalalignment='center', transform=ax1.transAxes, fontsize = 12,\
#         color = 'white')
#ax2.text(0.01, 0.1, 'Anisotropic A-D-C-B-A', horizontalalignment='left', \
#         verticalalignment='center', transform=ax2.transAxes, fontsize = 12,\
#         color = 'white')
#ax3.text(0.01, 0.1, 'Isotropic', horizontalalignment='left', \
#         verticalalignment='center', transform=ax3.transAxes, fontsize = 12,\
#         color = 'white')
ax1.set_title('Anisotropic A-B-C-D-A')
ax2.set_title('Anisotropic A-D-C-B-A')
ax3.set_title('Isotropic')
#ax1.legend()
     
for ax in (ax1,ax2,ax3):
    cc = ax.scatter(env.hexXmap, env.hexYmap, c = env.hexElevationMap, cmap = cm.gist_earth,s=60)
    ax.scatter(posA[0], posA[1], facecolor = 'r', edgecolor='black', s=60)
    ax.text(posA[0]-1, posA[1],'A',color = 'white')
    ax.scatter(posB[0], posB[1], facecolor = 'r', edgecolor='black', s=60)
    ax.text(posB[0]+1, posB[1],'B',color = 'white')
    ax.scatter(posC[0], posC[1], facecolor = 'r', edgecolor='black', s=60)
    ax.text(posC[0]-1, posC[1],'C',color = 'white')
    ax.scatter(posD[0], posD[1], facecolor = 'r', edgecolor='black', s=60)
    ax.text(posD[0]+1, posD[1],'D',color = 'white')
    ax.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
    ax.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
    ax.set_aspect('equal')
#cc = ax4.scatter(env.hexXmap, env.hexYmap, c = env.hexElevationMap, cmap = cm.gist_earth,s=60)
#ax4.plot(-1,-1,'lime')
#ax4.plot(-1,-1,'c')
#ax4.plot(-1,-1,'r')
#ax4.legend(('A-B-C-D-A', 'A-D-C-B-A', 'Isotropic equivalent'))
#plt.axis('off')
#cc.set_visible(False)
#fig.colorbar(cc, ax=ax4, loc='left')
#plt.subplots_adjust(left = 0.00, right = 0.1, bottom = 0.0, top = 0.1, wspace = -1.0, hspace = -1.0)
fig.tight_layout()

cbar_ax = fig.add_axes([0.05, 0.18, 0.2, 0.02])
legend_ax = fig.add_axes([0.55, 0.17, 0.4, 0.05])
legend_ax.plot(0,0,'r',label='$W_ϕ = 1$\n$W_θ = 1$')
legend_ax.plot(0,0,'m',label='$W_ϕ = 1 + 3 tan_{α}$\n$W_θ = 1$')
legend_ax.plot(0,0,'b',label='$W_ϕ = 1 + 6 tan_{α}$\n$W_θ = 1$')
legend_ax.plot(0,0,'c',label='$W_ϕ = 1 + 3 tan_{α}$\n$W_θ = 1 + 3 tan_{α}$')
legend_ax.plot(0,0,'g',label='$W_ϕ = 1 + 6 tan_{α}$\n$W_θ = 1 + 6 tan_{α}$')
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(b=False)
plt.axis('off')
cbar = fig.colorbar(cc, cax=cbar_ax, orientation = 'horizontal')
cbar.set_label('Elevation [m]',color='w')
cbar.ax.tick_params(colors='w')
cbar.outline.set_edgecolor('w')
#cbar.ax.yaxis.set_tick_params(color='w')
legend_ax.legend(ncol=5)
plt.show()


############### SHOW ROLL PATH ###############
#plt.style.use('default')
#plt.rcParams["font.family"] = "Constantia"
#plt.rcParams['mathtext.fontset'] = 'cm'
#plt.rcParams['mathtext.rm'] = 'serif'
#fig, axes = plt.subplots(figsize=(8, 8), \
#      nrows = 1, ncols = 3, \
#      sharex = 'all', sharey = 'all')
#ax1 = axes[0]
#ax2 = axes[1]
#ax3 = axes[2]
#
#
#plt.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.5, top = 1.0, wspace = 0.01, hspace = 0.05)
#fig.text(0.5, 0.005, 'X-axis [m]', ha='center')
#fig.text(0.005, 0.5, 'Y-axis [m]', va='center', rotation='vertical')
#
#showAnisoPath(env_CUAD02_scene01, 'c', ax1, ax2,'solid')
#showIsoPath(env_isoCUAD02_scene01, 'c', ax3, ax3,'solid')
#ax1.text(0.5, 0.95, 'Go traverses \n (CAMIS models)', horizontalalignment='center', \
#         verticalalignment='center', transform=ax1.transAxes, fontsize = 12,\
#         color = 'white')
#ax2.text(0.5, 0.95, 'Return traverses \n (CAMIS models)', horizontalalignment='center', \
#         verticalalignment='center', transform=ax2.transAxes, fontsize = 12,\
#         color = 'white')
#ax3.text(0.5, 0.95, 'Go and Return traverses \n (isotropic equivalent models)', horizontalalignment='center', \
#         verticalalignment='center', transform=ax3.transAxes, fontsize = 12,\
#         color = 'white')
##ax1.legend()
#     
#for ax in axes:
#    cc = ax.scatter(env.hexXmap, env.hexYmap, c = env.hexElevationMap, cmap = cm.gist_earth,s=20)
#    ax.scatter(posA[0], posA[1], facecolor = 'r', edgecolor='black', s=60)
#    ax.scatter(posB[0], posB[1], facecolor = 'r', edgecolor='black', s=60)
#    ax.scatter(posC[0], posC[1], facecolor = 'r', edgecolor='black', s=60)
#    ax.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
#    ax.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
#    ax.set_aspect('equal')
#fig.tight_layout()



############### ORIENTATION ###############

def showOrientation(env,ax,color,name):
    for i in range(8):
        if i == 0:
            offSet = 0
        else:
            offSet = offSet + env[i-1].pathTravDist[-1]
        if i == 0:
            pointName = 'Xa'
        elif i == 1 or i == 5:
            pointName = 'Xd'
        elif i == 3:
            pointName = 'Xb'
        elif i == 2 or i == 4:
            pointName = 'Xc'
        ax.axvline(x=offSet, color='w', linewidth=4)
        ax.axvline(x=offSet, color='steelblue', linewidth=1)
        ax.annotate(pointName,
                    xy=(offSet, 15),
                    xytext=(2, 0),
                    textcoords="offset points",
                    ha='left', va='bottom',color='steelblue')
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
    ax.axvspan(offSet + env[-1].pathTravDist[-1], 180, alpha=1.0, color='w')
    ax.axvline(x=offSet + env[-1].pathTravDist[-1],color='w', linewidth=4)
    ax.axvline(x=offSet + env[-1].pathTravDist[-1],color='steelblue', linewidth=1)
    ax.annotate('Xa',
                xy=(offSet + env[-1].pathTravDist[-1], 15),
                xytext=(2, 0),
                textcoords="offset points",
                ha='left', va='bottom',color='steelblue')
    ax.annotate(name,
                xy=(179, 5),
                xytext=(15, -15),  # 3 points vertical offset
                textcoords="offset points",
                ha='right', va='bottom',color=color)

plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, rowaxes = plt.subplots(figsize=(7, 7), nrows = 12, ncols = 1, \
     sharex = 'all', sharey = 'all')

rowaxes[0].set_visible(False)
rowaxes[1].set_title('Anisotropic')
rowaxes[7].set_title('Isotropic')
rowaxes[6].set_visible(False)
rowaxes[11].set_xlabel('Traversed Distance [m]')

plt.subplots_adjust(left = 0.085, right = 0.95, bottom = 0.075, top = 0.9, wspace = 0.1, hspace = 0.075)
fig.text(0.02, 0.5, 'Orientation Angle [deg]', va='center', rotation='vertical')
rowaxes[0].set_xlim([0,179])

showOrientation(env_CUAD01_scene01, rowaxes[1], 'r', '$W_ϕ = 1$\n$W_θ = 1$')
showOrientation(env_CUAD02_scene01, rowaxes[2], 'm', '$W_ϕ = 1 + 3 tan_{α}$\n$W_θ = 1$')
showOrientation(env_CUAD03_scene01, rowaxes[3], 'b', '$W_ϕ = 1 + 6 tan_{α}$\n$W_θ = 1$')
showOrientation(env_CUAD04_scene01, rowaxes[4], 'c', '$W_ϕ = 1 + 6 tan_{α}$\n$W_θ = 1 + 3 tan_{α}$')
showOrientation(env_CUAD05_scene01, rowaxes[5], 'g', '$W_ϕ = 1 + 6 tan_{α}$\n$W_θ = 1 + 6 tan_{α}$')
#showOrientation(env_CUAD06_scene01, rowaxes[5], 'c', 'CORI02')

showOrientation(env_isoCUAD01_scene01, rowaxes[7], 'r', '$W_ϕ = 1$\n$W_θ = 1$')
showOrientation(env_isoCUAD02_scene01, rowaxes[8], 'm', '$W_ϕ = 1 + 3 tan_{α}$\n$W_θ = 1$')
showOrientation(env_isoCUAD03_scene01, rowaxes[9], 'b', '$W_ϕ = 1 + 6 tan_{α}$\n$W_θ = 1$')
showOrientation(env_isoCUAD04_scene01, rowaxes[10], 'c', '$W_ϕ = 1 + 6 tan_{α}$\n$W_θ = 1 + 3 tan_{α}$')
showOrientation(env_isoCUAD05_scene01, rowaxes[11], 'g', '$W_ϕ = 1 + 6 tan_{α}$\n$W_θ = 1 + 6 tan_{α}$')
#showOrientation(env_isoCUAD06_scene01, rowaxes[11], 'y', 'isoCORI02')

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



## SHOW 3D PATHS

fig = plt.figure(figsize=plt.figaspect(0.5)) #Adjusts the aspect ratio and enlarges the figure (text does not enlarge)
ax = fig.gca(projection='3d')
XX = np.zeros_like(localDEM)
YY = np.zeros_like(localDEM)
for i in range(250):
    for j in range(250):
        XX[j,i] = offset[0] + i*0.1
        YY[j,i] = offset[1] + j*0.1

surf = ax.plot_surface(XX,YY,localDEM, cmap=cm.gist_earth,
                       linewidth=0, antialiased=False)
for i,anisomap in enumerate(env_CUAD02_scene01):
    if i==0 or i==2:
        ax.plot(anisomap.path[:,0]+offset[0],
                anisomap.path[:,1]+offset[1],
                anisomap.pathElevation,linestyle='solid')
    else:
        ax.plot(anisomap.path[:,0]+offset[0],
                anisomap.path[:,1]+offset[1],
                anisomap.pathElevation,linestyle='dashed')    
for i,anisomap in enumerate(env_isoCUAD02_scene01):
    if i==0 or i==2:
        ax.plot(anisomap.path[:,0]+offset[0],
                anisomap.path[:,1]+offset[1],
                anisomap.pathElevation,'r',linestyle='dotted')
    else:
        ax.plot(anisomap.path[:,0]+offset[0],
                anisomap.path[:,1]+offset[1],
                anisomap.pathElevation,'r',linestyle='dotted') 
MAX = 6
for direction in (-1, 1):
    for point in np.diag(direction * MAX * np.array([1,1,1])):
        ax.plot([point[0]+offset[0]+12.5], [point[1]+offset[1]+12.5], [point[2]+55.5], 'w')


################ ENERGY PLOT ##################
    
        
coeffsLabels = ['CUAD01','CUAD02','CUAD03','CUAD04','CUAD05']
anisoTotalCost = [0,0,0,0,0]
isoTotalCost = [0,0,0,0,0]
anisoTR = [0,0]
isoTR = [0,0]
for i in range(8):
    anisoTotalCost[0] = anisoTotalCost[0] + env_CUAD01_scene01[i].pathComputedTotalCost[-1]/3600.0
    anisoTotalCost[1] = anisoTotalCost[1] + env_CUAD02_scene01[i].pathComputedTotalCost[-1]/3600.0
    anisoTotalCost[2] = anisoTotalCost[2] + env_CUAD03_scene01[i].pathComputedTotalCost[-1]/3600.0
    anisoTotalCost[3] = anisoTotalCost[3] + env_CUAD04_scene01[i].pathComputedTotalCost[-1]/3600.0
    anisoTotalCost[4] = anisoTotalCost[4] + env_CUAD05_scene01[i].pathComputedTotalCost[-1]/3600.0
#    anisoTotalCost[5] = anisoTotalCost[5] + env_CUAD06_scene01[i].pathComputedTotalCost[-1]/3600.0
#    anisoTotalCost[6] = anisoTotalCost[6] + env_CUAD06_scene01[i].pathComputedTotalCost[-1]/3600.0
    isoTotalCost[0] = isoTotalCost[0] + env_isoCUAD01_scene01[i].pathComputedTotalCost[-1]/3600.0
    isoTotalCost[1] = isoTotalCost[1] + env_isoCUAD02_scene01[i].pathComputedTotalCost[-1]/3600.0
    isoTotalCost[2] = isoTotalCost[2] + env_isoCUAD03_scene01[i].pathComputedTotalCost[-1]/3600.0
    isoTotalCost[3] = isoTotalCost[3] + env_isoCUAD04_scene01[i].pathComputedTotalCost[-1]/3600.0
    isoTotalCost[4] = isoTotalCost[4] + env_isoCUAD05_scene01[i].pathComputedTotalCost[-1]/3600.0
#    isoTotalCost[5] = isoTotalCost[5] + env_isoCUAD06_scene01[i].pathComputedTotalCost[-1]/3600.0
#    isoTotalCost[6] = isoTotalCost[6] + env_isoCUAD06_scene01[i].pathComputedTotalCost[-1]/3600.0
    anisoTR[0] = anisoTR[0] + env_CUAD02_scene01[i].pathComputedTotalCostwithRisk[-1]/3600.0
    anisoTR[1] = anisoTR[1] + env_CUAD03_scene01[i].pathComputedTotalCostwithRisk[-1]/3600.0
    isoTR[0] = isoTR[0] + env_isoCUAD02_scene01[i].pathComputedTotalCostwithRisk[-1]/3600.0
    isoTR[1] = isoTR[1] + env_isoCUAD03_scene01[i].pathComputedTotalCostwithRisk[-1]/3600.0

x = np.arange(len(coeffsLabels))  # the label locations
x2 = np.arange(2)+1
width = 0.4  # the width of the bars

fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
rects3 = ax.bar(x2 - 0.45/2, anisoTR, 0.45, label='Isotropic (ρ = 0.8)', color='lime')
rects4 = ax.bar(x2 + 0.45/2, isoTR, 0.45, label='Isotropic (ρ = 0.8)', color='g')
rects1 = ax.bar(x - width/2, anisoTotalCost, width, label='Isotropic (ρ = 0.8)', color='r')
rects2 = ax.bar(x + width/2, isoTotalCost, width, label='Isotropic (ρ = 0.8)', color = 'b')

for i,rect in enumerate(rects1):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T) + '\nAh',
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, -3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='top',color = 'w')
for i,rect in enumerate(rects2):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T) + '\nAh',
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, -3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='top',color = 'w')
for i,rect in enumerate(rects2):
        isoT = rect.get_height()
        anisoT = rects1[i].get_height()
        gain = (isoT - anisoT)/isoT * 100
        ax.annotate('Gain = ' + '{0:.2f}'.format(gain) + '%',
                    xy=(rect.get_x(), isoT),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
for i,rect in enumerate(rects3):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T) + '\nAh',
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom',color = 'k')
for i,rect in enumerate(rects4):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T) + '\nAh',
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom',color = 'k')
for i,rect in enumerate(rects4):
        isoT = rect.get_height()
        anisoT = rects3[i].get_height()
        gain = (isoT - anisoT)/isoT * 100
        ax.annotate('Gain = ' + '{0:.2f}'.format(gain) + '%',
                    xy=(rect.get_x(), isoT),
                    xytext=(0, 25),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')   
        
ax.grid(True, which='both')   
#autolabel(rects1)
#autolabel(rects2)
ax.set_ylabel('Total Cost [Ah]')
ax.set_xlabel('CAMIS')
ax.set_xticks(x)
ax.set_xticklabels(coeffsLabels)
ax.legend(('Anisotropic Total Cost + Risk','Isotropic Total Cost + Risk', 'Anisotropic Total Cost','Isotropic Total Cost'))
plt.minorticks_on()  
plt.show()







fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
slopeArray01 = []
slopeArray03 = []
slopeArrayiso07 = []
for k,anisomap in enumerate(env_CUAD01_scene01):
    slopeArray01 = np.concatenate((slopeArray01,np.abs(180.0/np.pi*np.asarray(env_CUAD01_scene01[k].pathRoll))))
for k,anisomap in enumerate(env_CUAD03_scene01):
    slopeArray03 = np.concatenate((slopeArray03,np.abs(180.0/np.pi*np.asarray(env_CUAD03_scene01[k].pathRoll))))
for k,anisomap in enumerate(env_isoCUAD07_scene01):
    slopeArrayiso07 = np.concatenate((slopeArrayiso07,np.abs(180.0/np.pi*np.asarray(env_isoCUAD07_scene01[k].pathRoll))))
plt.hist([slopeArray03,slopeArrayiso07], bins = 100)





################ ORIENTATION PER TRAVERSE ANALYSIS ##################
def computeaboveSlope(env_input, linearGradient):
    aboveSlope = np.zeros_like(linearGradient)
    for i,slopelim in enumerate(linearGradient):
        for k,anisomap in enumerate(env_input):
            segmentArray = env_input[k].pathSegment
            slopeArray = np.abs(180.0/np.pi*np.asarray(env_input[k].pathSlope))
            for j, slope in enumerate(slopeArray):
                if (slope >= slopelim):
                    aboveSlope[i] = aboveSlope[i] + segmentArray[j]
    return aboveSlope
def computeaboveRoll(env_input, linearGradient):
    aboveRoll = np.zeros_like(linearGradient)
    for i,rolllim in enumerate(linearGradient):
        for k,anisomap in enumerate(env_input):
            segmentArray = env_input[k].pathSegment
            rollArray = np.abs(180.0/np.pi*np.asarray(env_input[k].pathRoll))
            for j, roll in enumerate(rollArray):
                if (roll >= rolllim):
                    aboveRoll[i] = aboveRoll[i] + segmentArray[j]
    return aboveRoll
def computeabovePositivePitch(env_input, linearGradient): #Descending
    abovePitch = np.zeros_like(linearGradient)
    for i,pitchlim in enumerate(linearGradient):
        for k,anisomap in enumerate(env_input):
            segmentArray = env_input[k].pathSegment
            pitchArray = 180.0/np.pi*np.asarray(env_input[k].pathPitch)
            for j, pitch in enumerate(pitchArray):
                if (pitch >= pitchlim):
                    abovePitch[i] = abovePitch[i] + segmentArray[j]
    return abovePitch
def computeunderNegativePitch(env_input, linearGradient): #Descending
    underPitch = np.zeros_like(linearGradient)
    for i,pitchlim in enumerate(linearGradient):
        for k,anisomap in enumerate(env_input):
            segmentArray = env_input[k].pathSegment
            pitchArray = 180.0/np.pi*np.asarray(env_input[k].pathPitch)
            for j, pitch in enumerate(pitchArray):
                if (pitch <= - pitchlim):
                    underPitch[i] = underPitch[i] + segmentArray[j]
    return underPitch
linearGradient = np.linspace(0,30,301) 
## ROLL ##  
fig, ax = plt.subplots(figsize=(4,4), constrained_layout=True)
ax.plot(linearGradient,computeaboveRoll(env_CUAD01_scene01, linearGradient))
ax.plot(linearGradient,computeaboveRoll(env_CUAD02_scene01, linearGradient))
ax.plot(linearGradient,computeaboveRoll(env_CUAD03_scene01, linearGradient))
ax.plot(linearGradient,computeaboveRoll(env_CUAD04_scene01, linearGradient))
ax.plot(linearGradient,computeaboveRoll(env_CUAD05_scene01, linearGradient))
ax.plot(linearGradient,computeaboveRoll(env_isoCUAD01_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeaboveRoll(env_isoCUAD02_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeaboveRoll(env_isoCUAD03_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeaboveRoll(env_isoCUAD04_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeaboveRoll(env_isoCUAD05_scene01, linearGradient), linestyle='dashed')
ax.legend(('CUAD01','CUAD02','CUAD03','CUAD04','CUAD05','CUAD06','CUAD07', \
           'isoCUAD01','isoCUAD02','isoCUAD03','isoCUAD04','isoCUAD05'))
ax.set_xlabel('Absolute Roll Threshold [degrees]')
ax.set_ylabel('Traversed distance above absolute roll threshold [m]')
ax.set_xlim([0,16])
## SLOPE ##
fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(linearGradient,computeaboveSlope(env_CUAD01_scene01, linearGradient))
ax.plot(linearGradient,computeaboveSlope(env_CUAD02_scene01, linearGradient))
ax.plot(linearGradient,computeaboveSlope(env_CUAD03_scene01, linearGradient))
ax.plot(linearGradient,computeaboveSlope(env_CUAD04_scene01, linearGradient))
ax.plot(linearGradient,computeaboveSlope(env_CUAD05_scene01, linearGradient))
ax.plot(linearGradient,computeaboveSlope(env_isoCUAD01_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeaboveSlope(env_isoCUAD02_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeaboveSlope(env_isoCUAD03_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeaboveSlope(env_isoCUAD04_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeaboveSlope(env_isoCUAD05_scene01, linearGradient), linestyle='dashed')
ax.legend(('CUAD01','CUAD02','CUAD03','CUAD04','CUAD05','CUAD06','CUAD07', \
           'isoCUAD01','isoCUAD02','isoCUAD03','isoCUAD04','isoCUAD05', \
           'isoCUAD06','isoCUAD07'))
ax.set_xlabel('Absolute Slope Threshold [degrees]')
ax.set_ylabel('Traversed distance above absolute slope threshold [m]')
ax.set_xlim([0,26])
## DESCENT PITCH ##
fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(linearGradient,computeabovePositivePitch(env_CUAD01_scene01, linearGradient))
ax.plot(linearGradient,computeabovePositivePitch(env_CUAD02_scene01, linearGradient))
ax.plot(linearGradient,computeabovePositivePitch(env_CUAD03_scene01, linearGradient))
ax.plot(linearGradient,computeabovePositivePitch(env_CUAD04_scene01, linearGradient))
ax.plot(linearGradient,computeabovePositivePitch(env_CUAD05_scene01, linearGradient))
ax.plot(linearGradient,computeabovePositivePitch(env_isoCUAD01_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeabovePositivePitch(env_isoCUAD02_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeabovePositivePitch(env_isoCUAD03_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeabovePositivePitch(env_isoCUAD04_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeabovePositivePitch(env_isoCUAD05_scene01, linearGradient), linestyle='dashed')
ax.legend(('CUAD01','CUAD02','CUAD03','CUAD04','CUAD05','CUAD06','CUAD07', \
           'isoCUAD01','isoCUAD02','isoCUAD03','isoCUAD04','isoCUAD05', \
           'isoCUAD06','isoCUAD07'))
ax.set_xlabel('Pitch Threshold [degrees]')
ax.set_ylabel('Traversed distance above pitch threshold [m]')
ax.set_xlim([0,26])
## ASCENT PITCH ##
fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(linearGradient,computeunderNegativePitch(env_CUAD01_scene01, linearGradient))
ax.plot(linearGradient,computeunderNegativePitch(env_CUAD02_scene01, linearGradient))
ax.plot(linearGradient,computeunderNegativePitch(env_CUAD03_scene01, linearGradient))
ax.plot(linearGradient,computeunderNegativePitch(env_CUAD04_scene01, linearGradient))
ax.plot(linearGradient,computeunderNegativePitch(env_CUAD05_scene01, linearGradient))
ax.plot(linearGradient,computeunderNegativePitch(env_isoCUAD01_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeunderNegativePitch(env_isoCUAD02_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeunderNegativePitch(env_isoCUAD03_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeunderNegativePitch(env_isoCUAD04_scene01, linearGradient), linestyle='dashed')
ax.plot(linearGradient,computeunderNegativePitch(env_isoCUAD05_scene01, linearGradient), linestyle='dashed')
ax.legend(('CUAD01','CUAD02','CUAD03','CUAD04','CUAD05','CUAD06','CUAD07', \
           'isoCUAD01','isoCUAD02','isoCUAD03','isoCUAD04','isoCUAD05', \
           'isoCUAD06','isoCUAD07'))
ax.set_xlabel('Pitch Threshold [degrees]')
ax.set_ylabel('Traversed distance under - pitch threshold [m]')
ax.set_xlim([0,26])



from pyproj import Proj, transform
f = open("test_CAMIS_iso_08_BC.gpx","w")
f.write('<?xml version="1.0" encoding="UTF-8" standalone="no" ?>\n')
f.write('<gpx ')
f.write('version="1.1"\n')
f.write('creator="UMA-SpaceRoboticsLab">\n')
f.write('<trk>\n')
f.write('<name>camisCuadrigaPath.gpx</name>\n')
f.write('<trkseg>\n')

D0 = 180.0/np.pi
l1 = offset[0]/D0
(l1*D0 + 183)/6

p = Proj(proj='utm', zone=30, ellps='WGS84')
for i,isomap in enumerate(env_isoCUAD07_scene01):
    if i == 1:
        currentPath = isomap.path
        for j in range(0,isomap.path.shape[0],20):
            pointX = isomap.path[j,0]+offset[0]
            pointY = isomap.path[j,1]+offset[1]           
#            lon,lat = transform(inProj,outProj,pointX,pointY)
            lon, lat = p(pointX, pointY, inverse=True)
            f.write('<trkpt lat="' + format(lat,'.12f') + '" lon="' + format(lon,'.12f') + '">\n')
            f.write(' <ele>0.0</ele>\n')
            f.write('  </trkpt>\n')
f.write('</trkseg>\n')
f.write('</trk>\n')
f.write('</gpx>\n')
f.close()
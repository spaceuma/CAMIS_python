# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 09:37:46 2020

@author: Richi
"""

import numpy as np

import data.cuadrigaData.cuadriga_reader as cr
from context import camis
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yaml
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import copy

import scipy.signal
from scipy.optimize import curve_fit

## =============================================================================
### IMPORT PLANNED PATHS
## =============================================================================
# Import 20201030execution.spydata using the variable explorer

def buildPathInfo(posX, posY, time, roll, pitch, current):
    pathInfo = np.concatenate(( posX - env.offset[0], \
                                posY - env.offset[1], \
                                time, roll, pitch, \
                                current )).reshape((6, len(posX)))
    for i in range(len(posX)-1,-1,-1):
        if (np.linalg.norm(pathInfo[0:2,i] - posA) < 2.8) or \
        (np.linalg.norm(pathInfo[0:2,i] - posB) < 2.8) or \
        (np.linalg.norm(pathInfo[0:2,i] - posC) < 2.8) or \
        (np.linalg.norm(pathInfo[0:2,i] - posD) < 2.8) or \
        (pathInfo[0,i] < posA[0]) or \
        (pathInfo[1,i] > posA[1]) or \
        (pathInfo[1,i] > posB[1]) or \
        (pathInfo[1,i] < 3.0):
            pathInfo = np.delete(pathInfo, i, axis=1)
    pathInfo[2] = pathInfo[2] - pathInfo[2][0]
    return pathInfo


### Anisotropic Model A2D results
posX_A2D, posY_A2D, _, Roll_A2D, Pitch_A2D, _, Current_A2D, _, Distance_A2D, Segment_A2D, _, \
 dTime_A2D, Time_A2D = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_06_43.txt',2)
posX_A2D_02, posY_A2D_02, _, Roll_A2D_02, Pitch_A2D_02, _, Current_A2D_02, _, Distance_A2D, Segment_A2D, _, \
 dTime_A2D, Time_A2D_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_19_51.txt',1)
path_aniso_A2D = buildPathInfo(posX_A2D, posY_A2D, Time_A2D, \
                               Roll_A2D, Pitch_A2D,Current_A2D)
path_aniso_A2D_02 = buildPathInfo(posX_A2D_02, posY_A2D_02, Time_A2D_02, \
                                  Roll_A2D_02, Pitch_A2D_02,Current_A2D_02)
fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_aniso_A2D[2], path_aniso_A2D[3])
ax.plot(path_aniso_A2D_02[2], path_aniso_A2D_02[3])


### Anisotropic Model D2A results
posX_D2A, posY_D2A, _, Roll_D2A, Pitch_D2A, _, Current_D2A, _, Distance_D2A, Segment_D2A, _, \
 dTime_D2A, Time_D2A = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_10_23.txt',1)
posX_D2A_02, posY_D2A_02, _, Roll_D2A_02, Pitch_D2A_02, _, Current_D2A_02, _, Distance_D2A, Segment_D2A, _, \
 dTime_D2A, Time_D2A_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_22_22.txt',1)
path_aniso_D2A = buildPathInfo(posX_D2A, posY_D2A, Time_D2A, \
                               Roll_D2A, Pitch_D2A,Current_D2A)
path_aniso_D2A_02 = buildPathInfo(posX_D2A_02, posY_D2A_02, Time_D2A_02, \
                                  Roll_D2A_02, Pitch_D2A_02, Current_D2A_02)
fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_aniso_D2A[2], path_aniso_D2A[3])
ax.plot(path_aniso_D2A_02[2], path_aniso_D2A_02[3])

### Isotropic Model A2D results
posX_isoA2D, posY_isoA2D, heading_A2D, Roll_isoA2D, Pitch_isoA2D, Yaw_A2D, \
Current_isoA2D, Speed_A2D, Distance_A2D, Segment_A2D, GPSspeed_A2D, \
 dTime_A2D, Time_isoA2D = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_13_37.txt',1)
posX_isoA2D_02, posY_isoA2D_02, heading_A2D, Roll_isoA2D_02, Pitch_isoA2D_02, Yaw_A2D, \
Current_isoA2D_02, Speed_A2D, Distance_A2D, Segment_A2D, GPSspeed_A2D, \
 dTime_A2D, Time_isoA2D_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_24_59.txt',1)
path_iso_A2D = buildPathInfo(posX_isoA2D, posY_isoA2D, Time_isoA2D, \
                             Roll_isoA2D, Pitch_isoA2D,Current_isoA2D)
path_iso_A2D_02 = buildPathInfo(posX_isoA2D_02, posY_isoA2D_02, Time_isoA2D_02, \
                             Roll_isoA2D_02, Pitch_isoA2D_02,Current_isoA2D_02)

fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_iso_A2D[2], path_iso_A2D[3])
ax.plot(path_iso_A2D_02[2], path_iso_A2D_02[3])

### Isotropic Model D2A results
posX_isoD2A, posY_isoD2A, heading, Roll_isoD2A, Pitch_isoD2A, Yaw, \
Current_isoD2A, Speed_D2A, Distance_D2A, Segment_D2A, GPSspeed, \
 dTime_D2A, Time_isoD2A = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_16_13.txt',1)
posX_isoD2A_02, posY_isoD2A_02, heading_02, Roll_isoD2A_02, Pitch_isoD2A_02, Yaw, \
Current_isoD2A_02, Speed_D2A, Distance_D2A, Segment_D2A, GPSspeed, \
 dTime_D2A, Time_isoD2A_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_27_26.txt',1)
path_iso_D2A = buildPathInfo(posX_isoD2A, posY_isoD2A, Time_isoD2A, \
                             Roll_isoD2A, Pitch_isoD2A, Current_isoD2A)
path_iso_D2A_02 = buildPathInfo(posX_isoD2A_02, posY_isoD2A_02, Time_isoD2A_02, \
                             Roll_isoD2A_02, Pitch_isoD2A_02,Current_isoD2A_02)

fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_iso_D2A[2], path_iso_D2A[3])
ax.plot(path_iso_D2A_02[2], path_iso_D2A_02[3])

### Anisotropic Model C2B results
posX_C2B, posY_C2B, heading, Roll_C2B, Pitch_C2B, Yaw, \
Current_C2B, Speed, Distance_C2B, Segment_C2B, GPSspeed, \
 dTime_C2B, Time_C2B = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_31_06.txt',1)
posX_C2B_02, posY_C2B_02, heading, Roll_C2B_02, Pitch_C2B_02, Yaw, \
Current_C2B_02, Speed, Distance_C2B, Segment_C2B, GPSspeed, \
 dTime_C2B, Time_C2B_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_41_07.txt',1)
path_aniso_C2B = buildPathInfo(posX_C2B, posY_C2B, Time_C2B, \
                               Roll_C2B, Pitch_C2B,Current_C2B)
path_aniso_C2B_02 = buildPathInfo(posX_C2B_02, posY_C2B_02, Time_C2B_02, \
                                  Roll_C2B_02, Pitch_C2B_02, Current_C2B_02)

### Anisotropic Model B2C results
posX_B2C, posY_B2C, heading, Roll_B2C, Pitch_B2C, Yaw, \
Current_B2C, Speed, Distance_B2C, Segment_B2C, GPSspeed, \
 dTime_B2C, Time_B2C = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_33_17.txt',1)
posX_B2C_02, posY_B2C_02, heading, Roll_B2C, Pitch_B2C, Yaw, \
Current_B2C, Speed, Distance_B2C, Segment_B2C, GPSspeed, \
 dTime_B2C, Time_B2C = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_43_28.txt',1)
 
### Isotropic Model C2B results
posX_isoC2B, posY_isoC2B, heading, Roll_isoC2B, Pitch_isoC2B, Yaw, \
Current_isoC2B, Speed, Distance_C2B, Segment_C2B, GPSspeed, \
 dTime_C2B, Time_isoC2B = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_35_30.txt',1)
posX_isoC2B_02, posY_isoC2B_02, heading, Roll_isoC2B_02, Pitch_isoC2B_02, Yaw, \
Current_isoC2B_02, Speed, Distance_C2B, Segment_C2B, GPSspeed, \
 dTime_C2B, Time_isoC2B_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_45_47.txt',1)
path_iso_C2B = buildPathInfo(posX_isoC2B, posY_isoC2B, Time_isoC2B, \
                               Roll_isoC2B, Pitch_isoC2B, Current_isoC2B)
path_iso_C2B_02 = buildPathInfo(posX_isoC2B_02, posY_isoC2B_02, Time_isoC2B_02, \
                                  Roll_isoC2B_02, Pitch_isoC2B_02, Current_isoC2B_02)
fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_aniso_C2B[2], path_aniso_C2B[3])
ax.plot(path_aniso_C2B_02[2], path_aniso_C2B_02[3])
ax.plot(path_iso_C2B[2], path_iso_C2B[3])
ax.plot(path_iso_C2B_02[2], path_iso_C2B_02[3])
ax.set_ylim([-15, 15])
ax.set_title('C2B')
 
### Isotropic Model B2C results
posX_isoB2C, posY_isoB2C, heading, Roll_B2C, Pitch_B2C, Yaw, \
Current_B2C, Speed, Distance_B2C, Segment_B2C, GPSspeed, \
 dTime_B2C, Time_B2C = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_38_48.txt',1)
posX_isoB2C_02, posY_isoB2C_02, heading, Roll_B2C, Pitch_B2C, Yaw, \
Current_B2C, Speed, Distance_B2C, Segment_B2C, GPSspeed, \
 dTime_B2C, Time_B2C = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_48_42.txt',1)
 
 
### Visualizing paths
env = copy.deepcopy(env_CUAD03_scene01[0])

plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, (ax1,ax2) = plt.subplots(figsize=(8, 4),nrows = 1, ncols = 2, constrained_layout=True)
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*env.hexSlopeMap, 
                 cmap="nipy_spectral",s=60.0)
cc = ax2.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*env.hexSlopeMap, 
                 cmap="nipy_spectral",s=60.0)
cbar = fig.colorbar(cc, ax=ax1,shrink=0.6)
cbar.set_label('Steepness (deg)')
ax1.set_xlim([ 0, 25 ])
ax1.set_ylim([ 0, 25 ])
ax2.set_xlim([ 0, 25 ])
ax2.set_ylim([ 0, 25 ])
ax1.scatter(path_aniso_A2D[0], path_aniso_A2D[1], facecolor = 'r',s=60)
ax1.scatter(path_aniso_A2D_02[0], path_aniso_A2D_02[1], facecolor = 'orange',s=60)
ax1.scatter(posX_D2A - env.offset[0], posY_D2A - env.offset[1], facecolor = 'b',s=60)
ax1.scatter(posX_D2A_02 - env.offset[0], posY_D2A_02 - env.offset[1], facecolor = 'c',s=60)
ax1.scatter(path_iso_A2D_02[0], path_iso_A2D_02[1], facecolor = 'lime',s=60)
ax1.scatter(posX_isoD2A - env.offset[0], posY_isoD2A - env.offset[1], facecolor = 'g',s=60)
ax1.scatter(posA[0], posA[1], facecolor = 'r', edgecolor='black', s=60)
circ1 = ax1.add_artist(plt.Circle((posA[0],posA[1]),2.8, color='w', alpha = 0.5))
ax1.scatter(posD[0], posD[1], facecolor = 'r', edgecolor='black', s=60)
circ2 = ax1.add_artist(plt.Circle((posD[0],posD[1]),2.8, color='w', alpha = 0.5))

ax2.scatter(path_aniso_C2B[0], path_aniso_C2B[1], facecolor = 'r',s=60)
ax2.scatter(path_aniso_C2B_02[0], path_aniso_C2B_02[1], facecolor = 'orange',s=60)
ax2.scatter(posX_B2C - env.offset[0], posY_B2C - env.offset[1], facecolor = 'b',s=60)
ax2.scatter(posX_B2C_02 - env.offset[0], posY_B2C_02 - env.offset[1], facecolor = 'c',s=60)
ax2.scatter(path_iso_C2B[0], path_iso_C2B[1], facecolor = 'lime',s=60)
ax2.scatter(posX_isoB2C - env.offset[0], posY_isoB2C - env.offset[1], facecolor = 'g',s=60)
ax2.scatter(posX_isoB2C_02 - env.offset[0], posY_isoB2C_02 - env.offset[1], facecolor = 'g',s=60)
ax2.scatter(posB[0], posB[1], facecolor = 'r', edgecolor='black', s=60)
circ3 = ax2.add_artist(plt.Circle((posB[0],posB[1]),2.8, color='w', alpha = 0.5))
ax2.scatter(posC[0], posC[1], facecolor = 'r', edgecolor='black', s=60)
circ4 = ax2.add_artist(plt.Circle((posC[0],posC[1]),2.8, color='w', alpha = 0.5))

def computeGradient(Roll, Pitch):
    deg2rad = np.pi/180
    rad2deg = 180/np.pi
    Gradient = np.zeros_like(Roll)
    for i, r in enumerate(Roll):
        Gradient[i] = rad2deg*np.arccos(np.cos(deg2rad*Roll[i])*np.cos(deg2rad*Pitch[i]))
    return Gradient






#expDistance = []
#expSegment = []
#expRoll = []
#expPitch = []
#expSlope = []
#expCurrent = []
#expTime = []
#expDTime = []
#

#
#
# 
#posX_D2C, posY_D2C, heading, Roll_D2C, Pitch_D2C, Yaw, \
#Current_D2C, Speed, Distance_D2C, Segment_D2C, GPSspeed, \
# dTime_D2C, Time_D2C = cr.readCuadrigaData('experimental_results/2020_10_13/2020_10_13_14_34_47.txt',1)
#Slope_D2C = computeGradient(Roll_D2C,Pitch_D2C)
#
#posX_C2B, posY_C2B, heading, Roll_C2B, Pitch_C2B, Yaw, \
#Current_C2B, Speed, Distance_C2B, Segment_C2B, GPSspeed, \
# dTime_C2B, Time_C2B = cr.readCuadrigaData('experimental_results/2020_10_13/2020_10_13_14_37_22.txt',1)
#Slope_C2B = computeGradient(Roll_C2B,Pitch_C2B)
#
#posX_B2C, posY_B2C, heading, Roll_B2C, Pitch_B2C, Yaw, \
#Current_B2C, Speed, Distance_B2C, Segment_B2C, GPSspeed, \
# dTime_B2C, Time_B2C = cr.readCuadrigaData('experimental_results/2020_10_13/2020_10_13_14_37_22.txt',3)
#Slope_B2C = computeGradient(Roll_B2C,Pitch_B2C)
#
#posX_C2D, posY_C2D, heading, Roll_C2D, Pitch_C2D, Yaw, \
#Current_C2D, Speed, Distance_C2D, Segment_C2D, GPSspeed, \
# dTime_C2D, Time_C2D = cr.readCuadrigaData('experimental_results/2020_10_13/2020_10_13_14_37_22.txt',4)
#Slope_C2D = computeGradient(Roll_C2D,Pitch_C2D)
#
#posX_D2A, posY_D2A, heading, Roll_D2A, Pitch_D2A, Yaw, \
#Current_D2A, Speed_D2A, Distance_D2A, Segment_D2A, GPSspeed, \
# dTime_D2A, Time_D2A = cr.readCuadrigaData('experimental_results/2020_10_13/2020_10_13_14_46_09.txt',5)
#Slope_D2A = computeGradient(Roll_D2A,Pitch_D2A)
#
#
#expDistance.append(Distance_A2D)
#expDistance.append(Distance_D2C)
#expDistance.append(Distance_C2B)
#expDistance.append(Distance_B2C)
#expDistance.append(Distance_C2D)
#expDistance.append(Distance_D2A)
#
#expSegment.append(Segment_A2D)
#expSegment.append(Segment_D2C)
#expSegment.append(Segment_C2B)
#expSegment.append(Segment_B2C)
#expSegment.append(Segment_C2D)
#expSegment.append(Segment_D2A)
#
#expTime.append(Time_A2D)
#expTime.append(Time_D2C)
#expTime.append(Time_C2B)
#expTime.append(Time_B2C)
#expTime.append(Time_C2D)
#expTime.append(Time_D2A)
#
#expDTime.append(dTime_A2D)
#expDTime.append(dTime_D2C)
#expDTime.append(dTime_C2B)
#expDTime.append(dTime_B2C)
#expDTime.append(dTime_C2D)
#expDTime.append(dTime_D2A)
#
#expRoll.append(Roll_A2D)
#expRoll.append(Roll_D2C)
#expRoll.append(Roll_C2B)
#expRoll.append(Roll_B2C)
#expRoll.append(Roll_C2D)
#expRoll.append(Roll_D2A)
#
#expPitch.append(Pitch_A2D)
#expPitch.append(Pitch_D2C)
#expPitch.append(Pitch_C2B)
#expPitch.append(Pitch_B2C)
#expPitch.append(Pitch_C2D)
#expPitch.append(Pitch_D2A)
#
#expSlope.append(Slope_A2D)
#expSlope.append(Slope_D2C)
#expSlope.append(Slope_C2B)
#expSlope.append(Slope_B2C)
#expSlope.append(Slope_C2D)
#expSlope.append(Slope_D2A)
#
#expCurrent.append(Current_A2D)
#expCurrent.append(Current_D2C)
#expCurrent.append(Current_C2B)
#expCurrent.append(Current_B2C)
#expCurrent.append(Current_C2D)
#expCurrent.append(Current_D2A)
#
#def showAnisoPath(mapList, color, ax1, ax2, mode):
#    for i,anisomap in enumerate(mapList):
#        if i < 3:
#            anisomap.showPath(fig,ax1,color,mode)
#        else:
#            anisomap.showPath(fig,ax2,color,mode)
#def showIsoPath(mapList, color, ax1, ax2, mode):
#    for i,anisomap in enumerate(mapList):
#        if i < 3:
#            anisomap.showPath(fig,ax1,color,mode)
#        else:
#            anisomap.showPath(fig,ax2,color,mode)
#
#plt.style.use('default')
#plt.rcParams["font.family"] = "Constantia"
#plt.rcParams['mathtext.fontset'] = 'cm'
#plt.rcParams['mathtext.rm'] = 'serif'
#fig, axes = plt.subplots(figsize=(8, 8), \
#      nrows = 1, ncols = 2, \
#      sharex = 'all', sharey = 'all')
#ax1 = axes[0]
#ax2 = axes[1]
#
#
#plt.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.5, top = 1.0, wspace = 0.01, hspace = 0.05)
#fig.text(0.5, 0.005, 'X-axis [m]', ha='center')
#fig.text(0.005, 0.5, 'Y-axis [m]', va='center', rotation='vertical')
#
#showAnisoPath(env_CUAD02_scene01, 'c', ax1, ax2,'solid')
#showAnisoPath(env_CUAD02_scene01, 'c', ax1, ax2,'solid')
#
##ax1.plot(posX_A2D - env.offset[0], posY_A2D - env.offset[1], 'r', linestyle='solid')
##ax1.plot(posX_D2C - env.offset[0], posY_D2C - env.offset[1], 'r', linestyle='solid')
##ax1.plot(posX_C2B - env.offset[0], posY_C2B - env.offset[1], 'r', linestyle='solid')
##ax2.plot(posX_B2C - env.offset[0], posY_B2C - env.offset[1], 'r', linestyle='solid')
##ax2.plot(posX_C2D - env.offset[0], posY_C2D - env.offset[1], 'r', linestyle='solid')
##ax2.plot(posX_D2A - env.offset[0], posY_D2A - env.offset[1], 'r', linestyle='solid')
#ax1.scatter(posX_A2D - env.offset[0], posY_A2D - env.offset[1], facecolor = 'r',s=60)
#ax1.scatter(posX_D2C - env.offset[0], posY_D2C - env.offset[1], facecolor = 'r',s=60)
#ax1.scatter(posX_C2B - env.offset[0], posY_C2B - env.offset[1], facecolor = 'r',s=60)
#ax2.scatter(posX_B2C - env.offset[0], posY_B2C - env.offset[1], facecolor = 'r',s=60)
#ax2.scatter(posX_C2D - env.offset[0], posY_C2D - env.offset[1], facecolor = 'r',s=60)
#ax2.scatter(posX_D2A - env.offset[0], posY_D2A - env.offset[1], facecolor = 'r',s=60)
#
#ax1.text(0.5, 0.95, 'Go traverses \n (CAMIS models)', horizontalalignment='center', \
#         verticalalignment='center', transform=ax1.transAxes, fontsize = 12,\
#         color = 'white')
#ax2.text(0.5, 0.95, 'Return traverses \n (CAMIS models)', horizontalalignment='center', \
#         verticalalignment='center', transform=ax2.transAxes, fontsize = 12,\
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
#
#
#
################ ORIENTATION ###############
#def showExperimentalOrientation(ax,color,name):
#    for i in range(6):
#        if i == 0:
#            offSet = 0
#        else:
#            offSet = offSet + expTime[i-1][-1]
#        if i == 0:
#            pointName = 'Xa'
#        elif i == 1 or i == 5:
#            pointName = 'Xd'
#        elif i == 3:
#            pointName = 'Xb'
#        elif i == 2 or i == 4:
#            pointName = 'Xc'
#        ax.axvline(x=offSet, color='w', linewidth=4)
#        ax.axvline(x=offSet, color='steelblue', linewidth=1)
#        ax.annotate(pointName,
#                    xy=(offSet, 15),
#                    xytext=(2, 0),
#                    textcoords="offset points",
#                    ha='left', va='bottom',color='steelblue')
#        ax.fill_between(offSet + np.asarray(expTime[i]), 0, \
#                  np.asarray(expSlope[i]), facecolor = 'y',\
#                  alpha = 0.5)
#        ax.plot(offSet + np.asarray(expTime[i]), \
#                      expSlope[i], 'y')
#        ax.fill_between(offSet + np.asarray(expTime[i]), 0, \
#                  -expSlope[i], facecolor = 'y',\
#                  alpha = 0.5)
#        ax.plot(offSet + np.asarray(expTime[i]), \
#                      -expSlope[i], 'y')
#        ax.fill_between(offSet + np.asarray(expTime[i]), 0, \
#                  expPitch[i], facecolor = 'c',\
#                  alpha = 0.5)
#        ax.plot(offSet + np.asarray(expTime[i]), \
#                      expPitch[i], 'c')
#        ax.fill_between(offSet + np.asarray(expTime[i]), 0, \
#                  expRoll[i], facecolor = 'r',\
#                  alpha = 0.5)
#        ax.plot(offSet + np.asarray(expTime[i]), \
#                      expRoll[i], 'r')
#    ax.axvspan(offSet + expTime[-1][-1], 180, alpha=1.0, color='w')
#    ax.axvline(x=offSet + expTime[-1][-1],color='w', linewidth=4)
#    ax.axvline(x=offSet + expTime[-1][-1],color='steelblue', linewidth=1)
#    ax.annotate('Xa',
#                xy=(offSet + expTime[-1][-1], 15),
#                xytext=(2, 0),
#                textcoords="offset points",
#                ha='left', va='bottom',color='steelblue')
#    ax.annotate(name,
#                xy=(179, 5),
#                xytext=(-4, 0),  # 3 points vertical offset
#                textcoords="offset points",
#                ha='right', va='bottom',color=color)
#
#
#def showOrientation(env,ax,color,name):
#    for i in range(6):
#        if i == 0:
#            offSet = 0
#        else:
#            offSet = offSet + env[i-1].pathTravDist[-1]/0.5
#        if i == 0:
#            pointName = 'Xa'
#        elif i == 1 or i == 5:
#            pointName = 'Xd'
#        elif i == 3:
#            pointName = 'Xb'
#        elif i == 2 or i == 4:
#            pointName = 'Xc'
#        ax.axvline(x=offSet, color='w', linewidth=4)
#        ax.axvline(x=offSet, color='steelblue', linewidth=1)
#        ax.annotate(pointName,
#                    xy=(offSet, 15),
#                    xytext=(2, 0),
#                    textcoords="offset points",
#                    ha='left', va='bottom',color='steelblue')
#        ax.fill_between(offSet + np.asarray(env[i].pathTravDist)/0.5, 0, \
#                  180.0/np.pi*np.asarray(env[i].pathSlope), facecolor = 'y',\
#                  alpha = 0.5)
#        ax.plot(offSet + np.asarray(env[i].pathTravDist)/0.5, \
#                      180.0/np.pi*np.asarray(env[i].pathSlope), 'y')
#        ax.fill_between(offSet + np.asarray(env[i].pathTravDist)/0.5, 0, \
#                  -180.0/np.pi*np.asarray(env[i].pathSlope), facecolor = 'y',\
#                  alpha = 0.5)
#        ax.plot(offSet + np.asarray(env[i].pathTravDist)/0.5, \
#                      -180.0/np.pi*np.asarray(env[i].pathSlope), 'y')
#        ax.fill_between(offSet + np.asarray(env[i].pathTravDist)/0.5, 0, \
#                  180.0/np.pi*np.asarray(env[i].pathPitch), facecolor = 'c',\
#                  alpha = 0.5)
#        ax.plot(offSet + np.asarray(env[i].pathTravDist)/0.5, \
#                      180.0/np.pi*np.asarray(env[i].pathPitch), 'c')
#        ax.fill_between(offSet + np.asarray(env[i].pathTravDist)/0.5, 0, \
#                  180.0/np.pi*env[i].pathRoll*np.sign(env[i].pathBeta), facecolor = 'r',\
#                  alpha = 0.5)
#        ax.plot(offSet + np.asarray(env[i].pathTravDist)/0.5, \
#                      180.0/np.pi*env[i].pathRoll*np.sign(env[i].pathBeta), 'r')
#    ax.axvspan(offSet + env[-1].pathTravDist[-1]/0.5, 180, alpha=1.0, color='w')
#    ax.axvline(x=offSet + env[-1].pathTravDist[-1]/0.5,color='w', linewidth=4)
#    ax.axvline(x=offSet + env[-1].pathTravDist[-1]/0.5,color='steelblue', linewidth=1)
#    ax.annotate('Xa',
#                xy=(offSet + env[-1].pathTravDist[-1]/0.5, 15),
#                xytext=(2, 0),
#                textcoords="offset points",
#                ha='left', va='bottom',color='steelblue')
#    ax.annotate(name,
#                xy=(179, 5),
#                xytext=(-4, 0),  # 3 points vertical offset
#                textcoords="offset points",
#                ha='right', va='bottom',color=color)
#
#plt.style.use('seaborn-darkgrid')
#plt.rcParams["font.family"] = "Constantia"
#plt.rcParams['mathtext.fontset'] = 'cm'
#plt.rcParams['mathtext.rm'] = 'serif'
#fig, rowaxes = plt.subplots(figsize=(10, 7), nrows = 5, ncols = 1, \
#     sharex = 'all', sharey = 'all')
#
#rowaxes[4].set_xlabel('Traversed Distance [m]')
#
#plt.subplots_adjust(left = 0.085, right = 0.95, bottom = 0.075, top = 0.9, wspace = 0.1, hspace = 0.075)
#fig.text(0.02, 0.5, 'Orientation Angle [deg]', va='center', rotation='vertical')
#rowaxes[0].set_xlim([0,400])
#
#showOrientation(env_CUAD01_scene01, rowaxes[0], 'r', 'CUAD')
#showOrientation(env_CUAD02_scene01, rowaxes[1], 'c', 'CORI')
#showExperimentalOrientation(rowaxes[2], 'c', 'Experimental')
#
#showOrientation(env_isoCUAD02_scene01, rowaxes[3], 'c', 'isoCORI')
#showOrientation(env_isoCUAD03_scene01, rowaxes[4], 'y', 'isoHighCost')
#
#plt.grid(b=True, which = 'minor')
#legend_ax = fig.add_axes([0.2, 0.92, 0.6, 0.07])
#plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
#plt.grid(b=False)
#xArray = np.arange(0, 2*np.pi, 0.01)
#legend_ax.plot(xArray, np.sin(xArray),'y')
#legend_ax.plot(xArray, -np.sin(xArray),'y')
#legend_ax.fill_between(xArray, -np.sin(xArray), \
#                  np.sin(xArray), facecolor = 'y',\
#                  alpha = 0.5)
#legend_ax.plot(xArray+2*np.pi, np.sin(xArray),'c')
#legend_ax.fill_between(xArray+2*np.pi, 0, \
#                  np.sin(xArray), facecolor = 'c',\
#                  alpha = 0.5)
#legend_ax.plot(xArray+4*np.pi, np.sin(xArray),'r')
#legend_ax.fill_between(xArray+4*np.pi, 0, \
#                  np.sin(xArray), facecolor = 'r',\
#                  alpha = 0.5)
#legend_ax.text(np.pi/2, 0.0, 'Steepness', horizontalalignment='center', \
#         verticalalignment='center', fontsize = 10,\
#         color = 'white')
#legend_ax.text(3*np.pi/2, -0.5, '-α', horizontalalignment='center', \
#         verticalalignment='center', fontsize = 12,\
#         color = 'white')
#legend_ax.text(3*np.pi/2, 0.5, 'α', horizontalalignment='center', \
#         verticalalignment='center', fontsize = 12,\
#         color = 'white')
#legend_ax.text(5*np.pi/2, 0.5, 'Pitch', horizontalalignment='center', \
#         verticalalignment='center', fontsize = 10,\
#         color = 'white')
#legend_ax.text(7*np.pi/2, -0.5, 'θ', horizontalalignment='center', \
#         verticalalignment='center', fontsize = 10,\
#         color = 'white')
#legend_ax.text(9*np.pi/2, 0.5, 'Roll', horizontalalignment='center', \
#         verticalalignment='center', fontsize = 10,\
#         color = 'white')
#legend_ax.text(11*np.pi/2, -0.5, 'Φ', horizontalalignment='center', \
#         verticalalignment='center', fontsize = 10,\
#         color = 'white')
#legend_ax.set_facecolor('w')
#
#
### CURRENT CONSUMPTION ##
#plt.style.use('seaborn-darkgrid')
#plt.rcParams["font.family"] = "Constantia"
#plt.rcParams['mathtext.fontset'] = 'cm'
#plt.rcParams['mathtext.rm'] = 'serif'
#fig, axes = plt.subplots(figsize=(8, 8), \
#      nrows = 1, ncols = 1, \
#      sharex = 'all', sharey = 'all')
#
##axes.fill_between(Distance_D2A, 0, \
##                  np.divide(np.asarray(Current_D2A),np.asarray(Speed_D2A)), facecolor = 'c',\
##                  alpha = 0.5)
#
##pathSegment = np.zeros_like(env_CUAD02_scene01[-1].pathSegment)
##for index, waypoint in enumerate(env_CUAD02_scene01[-1].path):
##            if index == 0:
##                pathSegment[index] = 0.0
##            else:
##                A = env_CUAD02_scene01[-1].path[index] - env_CUAD02_scene01[-1].path[index-1]
##                pathSegment[index] = np.linalg.norm(A)               
##travDist = np.cumsum(pathSegment)
#
##pathSegment = np.zeros_like(posX_D2A)
##for index, waypoint in enumerate(posX_D2A):
##            if index == 0:
##                pathSegment[index] = 0.0
##            else:
##                A = [posX_D2A[index] - posX_D2A[index-1],posY_D2A[index] - posY_D2A[index-1]]
##                pathSegment[index] = np.linalg.norm(A)               
##travDist = np.cumsum(pathSegment)
#
#axes.fill_between(np.cumsum(env_CUAD02_scene01[-1].pathSegment), 0, \
#                  np.asarray(env_CUAD02_scene01[-1].pathCost)*0.5, facecolor = 'r',\
#                  alpha = 0.5)
#axes.fill_between(Distance_D2A, 0, \
#                  np.asarray(Current_D2A), facecolor = 'r',\
#                  alpha = 0.5)
#ax.set_title('Current Consumption')
#
#
#
#
#
#fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
#linearGradient = np.linspace(0,30,301)
#aboveSlope02 = np.zeros_like(linearGradient)
#aboveSlopeExp = np.zeros_like(linearGradient)
#
#for i,slopelim in enumerate(linearGradient):
#    for k,anisomap in enumerate(env_CUAD02_scene01):
#        segmentArray = env_CUAD02_scene01[k].pathSegment
#        slopeArray = np.abs(180.0/np.pi*np.asarray(env_CUAD02_scene01[k].pathSlope))
#        for j, slope in enumerate(slopeArray):
#            if (slope >= slopelim):
#                aboveSlope02[i] = aboveSlope02[i] + segmentArray[j]/0.5 # 0.5 = Cuadriga speed
#                
#for i,slopelim in enumerate(linearGradient):
#    for k,slopeA in enumerate(expSlope):
#        timeArray = expDTime[k]
#        slopeArray = np.abs(slopeA)
#        for j, slope in enumerate(slopeArray):
#            if (slope >= slopelim):
#                aboveSlopeExp[i] = aboveSlopeExp[i] + timeArray[j] 
#                
#ax.fill_between(linearGradient, 0, \
#                aboveSlope02, facecolor = 'g',\
#                alpha = 0.1)
#ax.fill_between(linearGradient, 0, \
#                aboveSlopeExp, facecolor = 'r',\
#                alpha = 0.1)
#
#fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
#slopeArray = []
#slopeArray02 = []
#for i,s in enumerate(expSlope):
#    slopeArray = np.concatenate((slopeArray,s))
#for k,anisomap in enumerate(env_CUAD02_scene01):
#    slopeArray02 = np.concatenate((slopeArray02,np.abs(180.0/np.pi*np.asarray(env_CUAD02_scene01[k].pathRoll))))
#plt.hist([slopeArray,slopeArray02], bins = 100)
#
#
### ROLL ##
#linearGradient = np.linspace(0,30,301)
#def computeaboveRoll(env_input, linearGradient):
#    aboveRoll = np.zeros_like(linearGradient)
#    for i,rolllim in enumerate(linearGradient):
#        for k,anisomap in enumerate(env_input):
#            segmentArray = env_input[k].pathSegment
#            rollArray = np.abs(180.0/np.pi*np.asarray(env_input[k].pathRoll))
#            for j, roll in enumerate(rollArray):
#                if (roll >= rolllim):
#                    aboveRoll[i] = aboveRoll[i] + segmentArray[j]
#    return aboveRoll
#aboveRollexp = np.zeros_like(linearGradient)
#for i,rolllim in enumerate(linearGradient):
#    for k, rollArray in enumerate(expRoll):
#        for j, roll in enumerate(rollArray):
#            if (abs(roll) >= rolllim):
#                aboveRollexp[i] = aboveRollexp[i] + expSegment[k][j]
#
#fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
#ax.plot(linearGradient,computeaboveRoll(env_CUAD02_scene01, linearGradient))
#ax.plot(linearGradient,aboveRollexp)
#
#
#plt.show()
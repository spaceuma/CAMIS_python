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
from scipy.spatial.transform import Rotation as R
from scipy import interpolate

## =============================================================================
### IMPORT PLANNED PATHS
## =============================================================================
# Import 20201030execution.spydata using the variable explorer
env = env_CUAD03_scene01[0]

def buildPathInfo(posX, posY, time, roll, pitch, yaw, current, devAngle, fixAngle, I2Cmatrix):
    logic = np.ones(len(time))
    posX = posX - env.offset[0]
    posY = posY - env.offset[1]
#    roll = roll + 180.0
#    yaw = yaw + 180.0
    for i,x in enumerate(posX):
        if (np.linalg.norm([posX[i],posY[i]] - posA) < 2.8) or \
        (np.linalg.norm([posX[i],posY[i]] - posB) < 2.8) or \
        (np.linalg.norm([posX[i],posY[i]] - posC) < 2.8) or \
        (np.linalg.norm([posX[i],posY[i]] - posD) < 2.8) or \
        (posX[i] < posA[0]) or \
        (posY[i] > posA[1]) or \
        (posY[i] > posB[1]) or \
        (posY[i] < 3.0):
            logic[i] = False
        else:
            logic[i] = True
            
    posX = posX[logic == 1.0]
    posY = posY[logic == 1.0]
    time = time[logic == 1.0]
    roll = roll[logic == 1.0]
    pitch = pitch[logic == 1.0]
    yaw = yaw[logic == 1.0]
    current = current[logic == 1.0]
    time = time - time[0]
     
    pathInfo = np.concatenate(( posX, \
                                posY, \
                                time, roll, pitch, yaw, \
                                current, time, time, time)).reshape((10, len(posX)))
#    pathInfo[2] = pathInfo[2] - pathInfo[2][0]
    
#    I2Cmatrix = R.from_euler('Z', fixAngle, degrees=True) * R.from_euler('ZYX', [yaw[refPosIndex], \
#                                     pitch[refPosIndex], \
#                                     roll[refPosIndex]], degrees=True)
    imuRot = I2Cmatrix.as_matrix()
    I2C_Z = [imuRot[0][2], imuRot[1][2], imuRot[2][2]]
    for index, waypoint in enumerate(pathInfo[0]):
        # Slope compensation
        imuRot = R.from_euler('Z', fixAngle, degrees=True) * R.from_euler('ZYX', [yaw[index], pitch[index], roll[index]], degrees=True) * R.from_euler('Z', devAngle, degrees=True)
        imuRot = imuRot.as_matrix()
        imuZ = [imuRot[0][2], imuRot[1][2], imuRot[2][2]]
        pathInfo[9][index] = np.arccos(np.sum(np.dot(I2C_Z,imuZ)))*180.0/3.1416
        # Roll, Pitch compensation
        I2R = R.from_euler('Z', fixAngle, degrees=True) * R.from_euler('ZYX', [yaw[index], pitch[index], roll[index]], degrees=True) * R.from_euler('Z', devAngle, degrees=True)
        C2R = I2Cmatrix.inv()*I2R
        pathInfo[5][index],pathInfo[4][index],pathInfo[3][index] = \
        C2R.as_euler('ZYX', degrees=True)
        # Segment length and heading angle
        if index == 0:
            A = [pathInfo[0][1] - pathInfo[0][0], pathInfo[1][1] - pathInfo[1][0]]
            pathInfo[7][index] = np.linalg.norm(A) * 0.5
            pathInfo[8][index] = np.arctan2(A[1],A[0]) * 180.0/3.1416
        elif index == pathInfo[0].size - 1:
            A = [pathInfo[0][index] - pathInfo[0][index-1], pathInfo[1][index] - pathInfo[1][index-1]]
            pathInfo[7][index] = np.linalg.norm(A) * 0.5
            pathInfo[8][index] = np.arctan2(A[1],A[0]) * 180.0/3.1416
        else:
            A1 = [pathInfo[0][index] - pathInfo[0][index-1], pathInfo[1][index] - pathInfo[1][index-1]]
            A2 = [pathInfo[0][index+1] - pathInfo[0][index], pathInfo[1][index+1] - pathInfo[1][index]]
            pathInfo[7][index] = np.linalg.norm(A1) * 0.5 + np.linalg.norm(A2) * 0.5
            pathInfo[8][index] = np.arctan2(A1[1]+A2[1],A1[0]+A1[0]) * 180.0/3.1416
    Ximu = []
    Yimu = []
    Zimu = []
    for i,x in enumerate(roll):
        rot = R.from_euler('Z', fixAngle, degrees=True) * R.from_euler('ZYX', [yaw[i], pitch[i], roll[i]], degrees=True) * R.from_euler('Z', devAngle, degrees=True)
        Ximu.append(rot.apply(np.array([1, 0, 0]),inverse = False))
        Yimu.append(rot.apply(np.array([0, 1, 0]),inverse = False))
        Zimu.append(rot.apply(np.array([0, 0, 1]),inverse = False))
    I2Cmatrix = I2Cmatrix.inv().as_matrix()
    print(I2Cmatrix)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #ax.quiver(np.zeros(len(Ximu)), np.zeros(len(Ximu)), np.zeros(len(Ximu)), \
    #          vx, vy, vz, arrow_length_ratio = 0.01,alpha=.2)
    ax.quiver(0.0, 0.0, 0.0, \
              0.0, 0.0, 1.0, arrow_length_ratio = 0.05, color='c',linewidth = 5.0)
    ax.quiver(0.0, 0.0, 0.0, \
              0.0, -1.0, 0.0, arrow_length_ratio = 0.05, color='b',linewidth = 5.0)
    colors = np.arange(len(Zimu))
    ax.scatter([x[0] for i,x in enumerate(Ximu)],\
               [x[1] for i,x in enumerate(Ximu)],\
               [x[2] for i,x in enumerate(Ximu)],c=colors, cmap = 'hsv')
#    ax.scatter([x[0] for i,x in enumerate(Yimu)],\
#               [x[1] for i,x in enumerate(Yimu)],\
#               [x[2] for i,x in enumerate(Yimu)],c=colors)
    ax.scatter([x[0] for i,x in enumerate(Zimu)],\
               [x[1] for i,x in enumerate(Zimu)],\
               [x[2] for i,x in enumerate(Zimu)],c=colors, cmap = 'hsv')
    ax.quiver(0.0, 0.0, 0.0, \
          I2Cmatrix[0][0], I2Cmatrix[0][1], I2Cmatrix[0][2], arrow_length_ratio = 0.05, color='k',linewidth = 5.0)
    ax.quiver(0.0, 0.0, 0.0, \
          I2Cmatrix[1][0], I2Cmatrix[1][1], I2Cmatrix[1][2], arrow_length_ratio = 0.05, color='k',linewidth = 5.0)
    ax.quiver(0.0, 0.0, 0.0, \
          I2Cmatrix[2][0], I2Cmatrix[2][1], I2Cmatrix[2][2], arrow_length_ratio = 0.05, color='k',linewidth = 5.0)
    ax.scatter(-0.5,-0.5,0.0,c='g',alpha = 0.0)
    ax.scatter(0.5,0.5,0.0,c='g',alpha = 0.0)
    ax.scatter(0.0,0.0,1.0,c='g')
    return pathInfo

I2Cmatrix = R.from_euler('ZYX', [40.775069257688706 - 90.0,\
                                     -12.1770199928157, \
                                     0.0], degrees=True) * R.from_euler('Z', -40.775069257688706 + 90.0, degrees=True)

### A2D results
posX_A2D, posY_A2D, _, Roll_A2D, Pitch_A2D, Yaw_A2D, \
Current_A2D, _, Distance_A2D, Segment_A2D, _, dTime_A2D, Time_A2D, \
Slope_A2D = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_06_43.txt',2)
 
posX_A2D_02, posY_A2D_02, _, Roll_A2D_02, Pitch_A2D_02, Yaw_A2D_02, Current_A2D_02, _, Distance_A2D, Segment_A2D, _, \
 dTime_A2D, Time_A2D_02, Slope_A2D_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_19_51.txt',1)
I2Cmatrix_A2D = R.from_euler('Y', -2.0, degrees=True)*I2Cmatrix
I2Cmatrix_A2D = R.from_euler('X', -2.0, degrees=True)*I2Cmatrix_A2D
path_aniso_A2D = buildPathInfo(posX_A2D, posY_A2D, Time_A2D, \
                               Roll_A2D, Pitch_A2D, Yaw_A2D, Current_A2D, -20.0, 10.0, I2Cmatrix_A2D)
path_aniso_A2D_02 = buildPathInfo(posX_A2D_02, posY_A2D_02, Time_A2D_02, \
                                  Roll_A2D_02, Pitch_A2D_02, Yaw_A2D_02,Current_A2D_02, -20.0, 10.0, I2Cmatrix_A2D)

posX_isoA2D, posY_isoA2D, heading_A2D, Roll_isoA2D, Pitch_isoA2D, Yaw_A2D, \
Current_isoA2D, Speed_A2D, Distance_A2D, Segment_A2D, GPSspeed_A2D, \
 dTime_isoA2D, Time_isoA2D, Slope_isoA2D = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_13_37.txt',1)
posX_isoA2D_02, posY_isoA2D_02, heading_A2D, Roll_isoA2D_02, Pitch_isoA2D_02, Yaw_A2D_02, \
Current_isoA2D_02, Speed_A2D, Distance_A2D, Segment_A2D, GPSspeed_A2D, \
 dTime_isoA2D, Time_isoA2D_02, Slope_isoA2D_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_24_59.txt',1)
I2Cmatrix_isoA2D = R.from_euler('Y', -2.0, degrees=True)*I2Cmatrix
I2Cmatrix_isoA2D = R.from_euler('X', -2.0, degrees=True)*I2Cmatrix_isoA2D 
path_iso_A2D = buildPathInfo(posX_isoA2D, posY_isoA2D, Time_isoA2D, \
                             Roll_isoA2D, Pitch_isoA2D, Yaw_A2D, Current_isoA2D, -20.0, 10.0, I2Cmatrix_isoA2D)
path_iso_A2D_02 = buildPathInfo(posX_isoA2D_02, posY_isoA2D_02, Time_isoA2D_02, \
                             Roll_isoA2D_02, Pitch_isoA2D_02, Yaw_A2D_02, Current_isoA2D_02, -20.0, 10.0, I2Cmatrix_isoA2D)

fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_aniso_A2D[2], path_aniso_A2D[3], 'b')
ax.plot(path_aniso_A2D_02[2], path_aniso_A2D_02[3], 'b')
f1 = interpolate.interp1d(path_aniso_A2D[2], path_aniso_A2D[3])
f2 = interpolate.interp1d(path_aniso_A2D_02[2], path_aniso_A2D_02[3])
timeVec = np.arange(0,np.min((path_aniso_A2D[2][-1],path_aniso_A2D_02[2][-1])), 0.02)
ax.fill_between(timeVec, f1(timeVec), f2(timeVec), facecolor = 'b', alpha = 0.5)
ax.plot(path_iso_A2D[2], path_iso_A2D[3], 'orange')
ax.plot(path_iso_A2D_02[2], path_iso_A2D_02[3], 'orange')
f1 = interpolate.interp1d(path_iso_A2D[2], path_iso_A2D[3])
f2 = interpolate.interp1d(path_iso_A2D_02[2], path_iso_A2D_02[3])
timeVec = np.arange(0,np.min((path_iso_A2D[2][-1],path_iso_A2D_02[2][-1])), 0.02)
ax.fill_between(timeVec, f1(timeVec), f2(timeVec), facecolor = 'orange', alpha = 0.5)
ax.plot(np.asarray(env_CUAD03_scene01[4].pathTravDist)  / 0.5, \
        np.asarray(env_CUAD03_scene01[4].pathRoll)*180.0/3.1416, 'c')
ax.plot(np.asarray(env_isoCUAD03_scene01[4].pathTravDist) / 0.5, \
        np.asarray(env_isoCUAD03_scene01[4].pathRoll)*180.0/3.1416, 'r')
ax.set_xlabel('Traverse Time [s]')
ax.set_ylabel('Roll [deg]')

fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(Yaw_A2D)
#ax.plot(Yaw_D2A)
#fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
#ax.scatter(path_aniso_A2D[5],path_aniso_A2D[8])
#ax.scatter(path_aniso_A2D_02[5],path_aniso_A2D_02[8])
#ax.set_xlabel('Yaw - IMU reading')
#ax.set_ylabel('Yaw - path tangent')
#
fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_aniso_A2D[9])
ax.plot(path_aniso_A2D_02[9])
ax.plot(path_iso_A2D[9])
ax.plot(path_iso_A2D_02[9])
ax.set_ylabel('Slope_A2D')

fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(Slope_A2D)

### D2A results
posX_D2A, posY_D2A, _, Roll_D2A, Pitch_D2A, Yaw_D2A, Current_D2A, _, Distance_D2A, Segment_D2A, _, \
 dTime_D2A, Time_D2A, Slope_D2A = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_10_23.txt',1)
posX_D2A_02, posY_D2A_02, _, Roll_D2A_02, Pitch_D2A_02, Yaw_D2A_02, Current_D2A_02, _, Distance_D2A, Segment_D2A, _, \
 dTime_D2A, Time_D2A_02, Slope_D2A_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_22_22.txt',1)
I2Cmatrix_D2A = R.from_euler('Y', -4.0, degrees=True)*I2Cmatrix
I2Cmatrix_D2A = R.from_euler('X', -2.0, degrees=True)*I2Cmatrix_D2A 
path_aniso_D2A = buildPathInfo(posX_D2A, posY_D2A, Time_D2A, \
                               Roll_D2A, Pitch_D2A, Yaw_D2A, Current_D2A, 105.0, 160.0,I2Cmatrix_D2A)
path_aniso_D2A_02 = buildPathInfo(posX_D2A_02, posY_D2A_02, Time_D2A_02, \
                                  Roll_D2A_02, Pitch_D2A_02, Yaw_D2A_02, Current_D2A_02, 105.0,160.0,I2Cmatrix_D2A)



### Isotropic Model D2A results
posX_isoD2A, posY_isoD2A, heading, Roll_isoD2A, Pitch_isoD2A, Yaw_isoD2A, \
Current_isoD2A, Speed_D2A, Distance_D2A, Segment_D2A, GPSspeed, \
 dTime_D2A, Time_isoD2A, Slope_isoD2A = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_27_26.txt',1)

#posX_isoD2A, posY_isoD2A, heading, Roll_isoD2A, Pitch_isoD2A, Yaw_isoD2A, \
#Current_isoD2A, Speed_D2A, Distance_D2A, Segment_D2A, GPSspeed, \
# dTime_D2A, Time_isoD2A, Slope_isoD2A = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_16_13.txt',1)
#path_iso_D2A = buildPathInfo(posX_isoD2A, posY_isoD2A, Time_isoD2A, \
#                             Roll_isoD2A, Pitch_isoD2A, Yaw_isoD2A, Current_isoD2A, 90.0, 10)
posX_isoD2A_02, posY_isoD2A_02, heading_02, Roll_isoD2A_02, Pitch_isoD2A_02, Yaw_isoD2A_02, \
Current_isoD2A_02, Speed_D2A, Distance_D2A, Segment_D2A, GPSspeed, \
 dTime_D2A, Time_isoD2A_02, Slope_isoD2A_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_12_06_16.txt',1)
I2Cmatrix_isoD2A = R.from_euler('Y', -8.0, degrees=True)*I2Cmatrix
I2Cmatrix_isoD2A = R.from_euler('X', 4.0, degrees=True)*I2Cmatrix_isoD2A 
path_iso_D2A_02 = buildPathInfo(posX_isoD2A_02, posY_isoD2A_02, Time_isoD2A_02, \
                             Roll_isoD2A_02, Pitch_isoD2A_02, Yaw_isoD2A_02, Current_isoD2A_02, 100.0, -160.0,I2Cmatrix_isoD2A)
path_iso_D2A = buildPathInfo(posX_isoD2A, posY_isoD2A, Time_isoD2A, \
                             Roll_isoD2A, Pitch_isoD2A, Yaw_isoD2A, Current_isoD2A, 100.0, -160.0,I2Cmatrix_isoD2A)


fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_aniso_D2A[2], path_aniso_D2A[3], 'b')
ax.plot(path_aniso_D2A_02[2], path_aniso_D2A_02[3], 'b')
f1 = interpolate.interp1d(path_aniso_D2A[2], path_aniso_D2A[3])
f2 = interpolate.interp1d(path_aniso_D2A_02[2], path_aniso_D2A_02[3])
timeVec = np.arange(0,np.min((path_aniso_D2A[2][-1],path_aniso_D2A_02[2][-1])), 0.02)
ax.fill_between(timeVec, f1(timeVec), f2(timeVec), facecolor = 'b', alpha = 0.5)
ax.plot(path_iso_D2A[2], path_iso_D2A[3], 'orange')
ax.plot(path_iso_D2A_02[2], path_iso_D2A_02[3], 'orange')
f1 = interpolate.interp1d(path_iso_D2A[2], path_iso_D2A[3])
f2 = interpolate.interp1d(path_iso_D2A_02[2], path_iso_D2A_02[3])
timeVec = np.arange(0,np.min((path_iso_D2A[2][-1],path_iso_D2A_02[2][-1])), 0.02)
ax.fill_between(timeVec, f1(timeVec), f2(timeVec), facecolor = 'orange', alpha = 0.5)
ax.plot(np.asarray(env_CUAD03_scene01[3].pathTravDist)  / 0.5, \
        np.asarray(env_CUAD03_scene01[3].pathRoll)*180.0/3.1416, 'c')
ax.plot(np.asarray(env_isoCUAD03_scene01[3].pathTravDist) / 0.5, \
        np.asarray(env_isoCUAD03_scene01[3].pathRoll)*180.0/3.1416, 'r')
ax.set_xlabel('Traverse Time [s]')
ax.set_ylabel('Roll [deg]')

#fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
#ax.scatter(path_aniso_D2A[5],path_aniso_D2A[8])
#ax.set_xlabel('Yaw - IMU reading')
#ax.set_ylabel('Yaw - path tangent')
#
fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_aniso_D2A[2], path_aniso_D2A[9])
ax.plot(path_aniso_D2A_02[2], path_aniso_D2A_02[9])
ax.plot(path_iso_D2A[2], path_iso_D2A[9])
ax.plot(path_iso_D2A_02[2], path_iso_D2A_02[9])
ax.set_ylabel('Slope_D2A')

### Anisotropic Model C2B results
posX_C2B, posY_C2B, heading, Roll_C2B, Pitch_C2B, Yaw_C2B, \
Current_C2B, Speed, Distance_C2B, Segment_C2B, GPSspeed, \
 dTime_C2B, Time_C2B, Slope_C2B = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_31_06.txt',1)
posX_C2B_02, posY_C2B_02, heading, Roll_C2B_02, Pitch_C2B_02, Yaw_C2B_02, \
Current_C2B_02, Speed, Distance_C2B, Segment_C2B, GPSspeed, \
 dTime_C2B, Time_C2B_02, Slope_C2B_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_41_07.txt',1)
I2Cmatrix_C2B = R.from_euler('Y', -4.0, degrees=True)*I2Cmatrix
I2Cmatrix_C2B = R.from_euler('X', -2.0, degrees=True)*I2Cmatrix_C2B 
path_aniso_C2B = buildPathInfo(posX_C2B, posY_C2B, Time_C2B, \
                               Roll_C2B, Pitch_C2B, Yaw_C2B,Current_C2B, 180.0, 90.0, I2Cmatrix_C2B)
path_aniso_C2B_02 = buildPathInfo(posX_C2B_02, posY_C2B_02, Time_C2B_02, \
                                  Roll_C2B_02, Pitch_C2B_02, Yaw_C2B_02, Current_C2B_02, 180.0, 90.0, I2Cmatrix_C2B)

posX_isoC2B, posY_isoC2B, heading, Roll_isoC2B, Pitch_isoC2B, Yaw_isoC2B, \
Current_isoC2B, Speed, Distance_C2B, Segment_C2B, GPSspeed, \
 dTime_C2B, Time_isoC2B, Slope_isoC2B = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_35_30.txt',1)
posX_isoC2B_02, posY_isoC2B_02, heading, Roll_isoC2B_02, Pitch_isoC2B_02, Yaw_isoC2B_02, \
Current_isoC2B_02, Speed, Distance_C2B, Segment_C2B, GPSspeed, \
 dTime_C2B, Time_isoC2B_02, Slope_isoC2B_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_45_47.txt',1)
I2Cmatrix_isoC2B = R.from_euler('Y', -4.0, degrees=True)*I2Cmatrix
I2Cmatrix_isoC2B = R.from_euler('X', 2.0, degrees=True)*I2Cmatrix_isoC2B 
path_iso_C2B = buildPathInfo(posX_isoC2B, posY_isoC2B, Time_isoC2B, \
                               Roll_isoC2B, Pitch_isoC2B, Yaw_isoC2B, Current_isoC2B, 30.0, 90.0, I2Cmatrix_isoC2B)
path_iso_C2B_02 = buildPathInfo(posX_isoC2B_02, posY_isoC2B_02, Time_isoC2B_02, \
                                  Roll_isoC2B_02, Pitch_isoC2B_02, Yaw_isoC2B_02, Current_isoC2B_02, 30.0, 90.0, I2Cmatrix_isoC2B)

fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_aniso_C2B[2], path_aniso_C2B[3], 'b')
ax.plot(path_aniso_C2B_02[2], path_aniso_C2B_02[3], 'b')
f1 = interpolate.interp1d(path_aniso_C2B[2], path_aniso_C2B[3])
f2 = interpolate.interp1d(path_aniso_C2B_02[2], path_aniso_C2B_02[3])
timeVec = np.arange(0,np.min((path_aniso_C2B[2][-1],path_aniso_C2B_02[2][-1])), 0.02)
ax.fill_between(timeVec, f1(timeVec), f2(timeVec), facecolor = 'b', alpha = 0.5)
ax.plot(path_iso_C2B[2], path_iso_C2B[3], 'orange')
ax.plot(path_iso_C2B_02[2], path_iso_C2B_02[3], 'orange')
f1 = interpolate.interp1d(path_iso_C2B[2], path_iso_C2B[3])
f2 = interpolate.interp1d(path_iso_C2B_02[2], path_iso_C2B_02[3])
timeVec = np.arange(0,np.min((path_iso_C2B[2][-1],path_iso_C2B_02[2][-1])), 0.02)
ax.fill_between(timeVec, f1(timeVec), f2(timeVec), facecolor = 'orange', alpha = 0.5)
ax.plot(np.asarray(env_CUAD03_scene01[6].pathTravDist)  / 0.5, \
        np.asarray(env_CUAD03_scene01[6].pathRoll)*180.0/3.1416, 'c')
ax.plot(np.asarray(env_isoCUAD03_scene01[6].pathTravDist) / 0.5, \
        np.asarray(env_isoCUAD03_scene01[6].pathRoll)*180.0/3.1416, 'r')
ax.set_xlabel('Traverse Time [s]')
ax.set_ylabel('Roll [deg]')

fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_aniso_C2B[9])
ax.plot(path_aniso_C2B_02[9])
ax.plot(path_iso_C2B[9])
ax.plot(path_iso_C2B_02[9])
ax.set_ylabel('Slope_C2B')


 
### Anisotropic Model B2C results
posX_B2C, posY_B2C, heading, Roll_B2C, Pitch_B2C, Yaw_B2C, \
Current_B2C, Speed, Distance_B2C, Segment_B2C, GPSspeed, \
 dTime_B2C, Time_B2C, Slope_B2C = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_33_17.txt',1)
posX_B2C_02, posY_B2C_02, heading, Roll_B2C_02, Pitch_B2C_02, Yaw_B2C_02, \
Current_B2C_02, Speed, Distance_B2C, Segment_B2C, GPSspeed, \
 dTime_B2C, Time_B2C_02, Slope_B2C_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_43_28.txt',1)
I2Cmatrix_B2C = R.from_euler('Y', -6.0, degrees=True)*I2Cmatrix
I2Cmatrix_B2C = R.from_euler('X', -1.0, degrees=True)*I2Cmatrix_B2C
path_aniso_B2C = buildPathInfo(posX_B2C, posY_B2C, Time_B2C, \
                               Roll_B2C, Pitch_B2C, Yaw_B2C, Current_B2C, -140.0, -20.0, I2Cmatrix_B2C)
path_aniso_B2C_02 = buildPathInfo(posX_B2C_02, posY_B2C_02, Time_B2C_02, \
                                  Roll_B2C_02, Pitch_B2C_02, Yaw_B2C_02, Current_B2C_02, -140.0, -20.0, I2Cmatrix_B2C)

### Isotropic Model B2C results
posX_isoB2C, posY_isoB2C, heading, Roll_isoB2C, Pitch_isoB2C, Yaw_isoB2C, \
Current_isoB2C, Speed, Distance_B2C, Segment_B2C, GPSspeed, \
 dTime_isoB2C, Time_isoB2C, Slope_isoB2C = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_38_48.txt',1)
posX_isoB2C_02, posY_isoB2C_02, heading, Roll_isoB2C_02, Pitch_isoB2C_02, Yaw_isoB2C_02, \
Current_isoB2C_02, Speed, Distance_B2C, Segment_B2C, GPSspeed, \
 dTime_B2C, Time_isoB2C_02, Slope_isoB2C_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_48_42.txt',1)
I2Cmatrix_isoB2C = R.from_euler('Y', -4.0, degrees=True)*I2Cmatrix
I2Cmatrix_isoB2C = R.from_euler('X', -2.0, degrees=True)*I2Cmatrix_isoB2C
path_iso_B2C = buildPathInfo(posX_isoB2C, posY_isoB2C, Time_isoB2C, \
                               Roll_isoB2C, Pitch_isoB2C, Yaw_isoB2C,Current_isoB2C, 40.0, -20.0, I2Cmatrix_isoB2C)
path_iso_B2C_02 = buildPathInfo(posX_isoB2C_02, posY_isoB2C_02, Time_isoB2C_02, \
                                  Roll_isoB2C_02, Pitch_isoB2C_02, Yaw_isoB2C_02, Current_isoB2C_02, 0.0, -20.0, I2Cmatrix_isoB2C) 
 
fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_aniso_B2C[2], path_aniso_B2C[3], 'b')
ax.plot(path_aniso_B2C_02[2], path_aniso_B2C_02[3], 'b')
f1 = interpolate.interp1d(path_aniso_B2C[2], path_aniso_B2C[3])
f2 = interpolate.interp1d(path_aniso_B2C_02[2], path_aniso_B2C_02[3])
timeVec = np.arange(0,np.min((path_aniso_B2C[2][-1],path_aniso_B2C_02[2][-1])), 0.02)
ax.fill_between(timeVec, f1(timeVec), f2(timeVec), facecolor = 'b', alpha = 0.5)
ax.plot(path_iso_B2C[2], path_iso_B2C[3], 'orange')
ax.plot(path_iso_B2C_02[2], path_iso_B2C_02[3], 'orange')
f1 = interpolate.interp1d(path_iso_B2C[2], path_iso_B2C[3])
f2 = interpolate.interp1d(path_iso_B2C_02[2], path_iso_B2C_02[3])
timeVec = np.arange(0,np.min((path_iso_B2C[2][-1],path_iso_B2C_02[2][-1])), 0.02)
ax.fill_between(timeVec, f1(timeVec), f2(timeVec), facecolor = 'orange', alpha = 0.5)
ax.plot(np.asarray(env_CUAD03_scene01[1].pathTravDist)  / 0.5, \
        np.asarray(env_CUAD03_scene01[1].pathRoll)*180.0/3.1416, 'c')
ax.plot(np.asarray(env_isoCUAD03_scene01[1].pathTravDist) / 0.5, \
        np.asarray(env_isoCUAD03_scene01[1].pathRoll)*180.0/3.1416, 'r')
ax.set_xlabel('Traverse Time [s]')
ax.set_ylabel('Roll [deg]')

fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.scatter(path_aniso_B2C[5],path_aniso_B2C[8])
ax.set_xlabel('Yaw - IMU reading')
ax.set_ylabel('Yaw - path tangent')

fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_aniso_B2C[9])
ax.plot(path_aniso_B2C_02[9])
ax.plot(path_iso_B2C[9])
ax.plot(path_iso_B2C_02[9])
ax.set_ylabel('Slope_B2C')

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
#ax1.scatter(path_aniso_A2D_03[0], path_aniso_A2D_03[1], facecolor = 'k',s=60)
ax1.scatter(path_aniso_D2A[0], path_aniso_D2A[1], facecolor = 'b',s=60)
ax1.scatter(posX_D2A_02 - env.offset[0], posY_D2A_02 - env.offset[1], facecolor = 'c',s=60)
ax1.scatter(path_iso_A2D_02[0], path_iso_A2D_02[1], facecolor = 'lime',s=60)
ax1.scatter(path_iso_D2A[0], path_iso_D2A[1], facecolor = 'g',s=60)
#ax1.scatter(path_iso_D2A_02[0], path_iso_D2A_02[1], facecolor = 'g',s=60)
ax1.scatter(posA[0], posA[1], facecolor = 'r', edgecolor='black', s=60)
circ1 = ax1.add_artist(plt.Circle((posA[0],posA[1]),2.8, color='w', alpha = 0.5))
ax1.scatter(posD[0], posD[1], facecolor = 'r', edgecolor='black', s=60)
circ2 = ax1.add_artist(plt.Circle((posD[0],posD[1]),2.8, color='w', alpha = 0.5))

ax2.scatter(path_aniso_C2B[0], path_aniso_C2B[1], facecolor = 'r',s=60)
ax2.scatter(path_aniso_C2B_02[0], path_aniso_C2B_02[1], facecolor = 'orange',s=60)
#ax2.scatter(posX_B2C - env.offset[0], posY_B2C - env.offset[1], facecolor = 'b',s=60)
ax2.scatter(path_aniso_B2C[0], path_aniso_B2C[1], facecolor = 'b',s=60)
ax2.scatter(path_aniso_B2C_02[0], path_aniso_B2C_02[1], facecolor = 'c',s=60)
ax2.scatter(path_iso_C2B[0], path_iso_C2B[1], facecolor = 'lime',s=60)
ax2.scatter(path_iso_C2B_02[0], path_iso_C2B_02[1], facecolor = 'lime',s=60)
ax2.scatter(path_iso_B2C[0], path_iso_B2C[1], facecolor = 'g',s=60)
ax2.scatter(path_iso_B2C_02[0], path_iso_B2C_02[1], facecolor = 'g',s=60)
ax2.scatter(posB[0], posB[1], facecolor = 'r', edgecolor='black', s=60)
circ3 = ax2.add_artist(plt.Circle((posB[0],posB[1]),2.8, color='w', alpha = 0.5))
ax2.scatter(posC[0], posC[1], facecolor = 'r', edgecolor='black', s=60)
circ4 = ax2.add_artist(plt.Circle((posC[0],posC[1]),2.8, color='w', alpha = 0.5))
ax1.plot(env_CUAD03_scene01[4].path[:,0], env_CUAD03_scene01[4].path[:,1],'k')
ax1.plot(env_CUAD03_scene01[3].path[:,0], env_CUAD03_scene01[3].path[:,1],'k')
ax2.plot(env_CUAD03_scene01[6].path[:,0], env_CUAD03_scene01[6].path[:,1],'k')
ax2.plot(env_CUAD03_scene01[1].path[:,0], env_CUAD03_scene01[1].path[:,1],'k')
ax1.plot(env_isoCUAD07_scene01[4].path[:,0], env_isoCUAD07_scene01[4].path[:,1],'k')
ax2.plot(env_isoCUAD07_scene01[6].path[:,0], env_isoCUAD07_scene01[6].path[:,1],'k')

def computeGradient(Roll, Pitch):
    deg2rad = np.pi/180
    rad2deg = 180/np.pi
    Gradient = np.zeros_like(Roll)
    for i, r in enumerate(Roll):
        Gradient[i] = rad2deg*np.arccos(np.cos(deg2rad*Roll[i])*np.cos(deg2rad*Pitch[i]))
    return Gradient


def computeaboveRoll(pathInfo, linearGradient):
    aboveRoll = np.zeros_like(linearGradient)
    for i,rolllim in enumerate(linearGradient):
        dTarray = np.append(0,np.diff(pathInfo[2]))
#        dTarray = pathInfo[7]
        rollArray = np.abs(pathInfo[3])
        for j, roll in enumerate(rollArray):
            if (roll >= rolllim):
                aboveRoll[i] = aboveRoll[i] + dTarray[j]
    return aboveRoll

def computeabovePitch(pathInfo, linearGradient):
    aboveRoll = np.zeros_like(linearGradient)
    for i,rolllim in enumerate(linearGradient):
        dTarray = np.append(0,np.diff(pathInfo[2]))
#        dTarray = pathInfo[7]
        rollArray = np.abs(pathInfo[4])
        for j, roll in enumerate(rollArray):
            if (roll >= rolllim):
                aboveRoll[i] = aboveRoll[i] + dTarray[j]
    return aboveRoll

linearGradient = np.linspace(0,30,301) 
## ROLL ##  
fig, ax = plt.subplots(figsize=(6,6), constrained_layout=True)
ax.plot(linearGradient,computeaboveRoll(path_aniso_A2D, linearGradient), 'b')
ax.plot(linearGradient,computeaboveRoll(path_aniso_A2D_02, linearGradient), 'c')
ax.plot(linearGradient,computeaboveRoll(path_iso_A2D, linearGradient),'r')
ax.plot(linearGradient,computeaboveRoll(path_iso_A2D_02, linearGradient),'orange')
ax.legend(('A2D - Anisotropic', 'A2D - Anisotropic', 'A2D - Isotropic', 'A2D - Isotropic'))
ax.set_xlabel('Absolute Roll Threshold [degrees]')
ax.set_ylabel('Distance above absolute roll threshold [m]')
ax.set_xlim([0,16])
## A2D
fig, ax = plt.subplots(figsize=(6,6), constrained_layout=True)
a1 = ax.fill_between(linearGradient, \
                  computeaboveRoll(path_aniso_A2D, linearGradient), \
                  computeaboveRoll(path_aniso_A2D_02, linearGradient), \
                  facecolor = 'b',\
                  alpha = 0.5)

a2 = ax.fill_between(linearGradient, \
                  computeaboveRoll(path_iso_A2D, linearGradient), \
                  computeaboveRoll(path_iso_A2D_02, linearGradient), \
                  facecolor = 'orange',\
                  alpha = 0.5)
pp = computeaboveRoll(path_aniso_A2D, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient,pp, 'b')
pp = computeaboveRoll(path_aniso_A2D_02, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient, pp, 'b')
pp = computeaboveRoll(path_iso_A2D, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient,pp, 'orange')
pp = computeaboveRoll(path_iso_A2D_02, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient, pp, 'orange')
ax.legend((a1,a2),('A2D - Anisotropic', 'A2D - Isotropic'))
ax.set_xlabel('Absolute Roll Threshold [degrees]')
ax.set_ylabel('Time above absolute roll threshold [m]')
ax.set_xlim([0,16])

## D2A
fig, ax = plt.subplots(figsize=(6,6), constrained_layout=True)
a1 = ax.fill_between(linearGradient, \
                  computeaboveRoll(path_aniso_D2A, linearGradient), \
                  computeaboveRoll(path_aniso_D2A_02, linearGradient), \
                  facecolor = 'b',\
                  alpha = 0.5)

a2 = ax.fill_between(linearGradient, \
                  computeaboveRoll(path_iso_D2A, linearGradient), \
                  computeaboveRoll(path_iso_D2A_02, linearGradient), \
                  facecolor = 'orange',\
                  alpha = 0.5)
pp = computeaboveRoll(path_aniso_D2A, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient,pp, 'b')
pp = computeaboveRoll(path_aniso_D2A_02, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient, pp, 'b')
pp = computeaboveRoll(path_iso_D2A, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient,pp, 'orange')
pp = computeaboveRoll(path_iso_D2A_02, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient, pp, 'orange')
ax.legend((a1,a2),('D2A - Anisotropic', 'D2A - Isotropic'))
ax.set_xlabel('Absolute Roll Threshold [degrees]')
ax.set_ylabel('Time above absolute roll threshold [m]')
ax.set_xlim([0,16])

## C2B
fig, ax = plt.subplots(figsize=(6,6), constrained_layout=True)
a1 = ax.fill_between(linearGradient, \
                  computeaboveRoll(path_aniso_C2B, linearGradient), \
                  computeaboveRoll(path_aniso_C2B_02, linearGradient), \
                  facecolor = 'b',\
                  alpha = 0.5)

a2 = ax.fill_between(linearGradient, \
                  computeaboveRoll(path_iso_C2B, linearGradient), \
                  computeaboveRoll(path_iso_C2B_02, linearGradient), \
                  facecolor = 'orange',\
                  alpha = 0.5)
pp = computeaboveRoll(path_aniso_C2B, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient,pp, 'b')
pp = computeaboveRoll(path_aniso_C2B_02, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient, pp, 'b')
pp = computeaboveRoll(path_iso_C2B, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient,pp, 'orange')
pp = computeaboveRoll(path_iso_C2B_02, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient, pp, 'orange')
ax.legend((a1,a2),('C2B - Anisotropic', 'C2B - Isotropic'))
ax.set_xlabel('Absolute Roll Threshold [degrees]')
ax.set_ylabel('Time above absolute roll threshold [m]')
ax.set_xlim([0,16])


#B2C
fig, ax = plt.subplots(figsize=(6,6), constrained_layout=True)
a1 = ax.fill_between(linearGradient, \
                  computeaboveRoll(path_aniso_B2C, linearGradient), \
                  computeaboveRoll(path_aniso_B2C_02, linearGradient), \
                  facecolor = 'b',\
                  alpha = 0.5)

a2 = ax.fill_between(linearGradient, \
                  computeaboveRoll(path_iso_B2C, linearGradient), \
                  computeaboveRoll(path_iso_B2C_02, linearGradient), \
                  facecolor = 'orange',\
                  alpha = 0.5)
pp = computeaboveRoll(path_aniso_B2C, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient,pp, 'b')
pp = computeaboveRoll(path_aniso_B2C_02, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient, pp, 'b')
pp = computeaboveRoll(path_iso_B2C, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient,pp, 'orange')
pp = computeaboveRoll(path_iso_B2C_02, linearGradient)
pp[np.where(pp == 0)] = np.nan
ax.plot(linearGradient, pp, 'orange')
ax.legend((a1,a2),('B2C - Anisotropic', 'B2C - Isotropic'))
ax.set_xlabel('Absolute Roll Threshold [degrees]')
ax.set_ylabel('Time above absolute roll threshold [m]')
ax.set_xlim([0,16])

#C2B2C
#fig, ax = plt.subplots(figsize=(6,6), constrained_layout=True)
#a1 = ax.fill_between(linearGradient, \
#                  computeaboveRoll(path_aniso_C2B, linearGradient) + computeaboveRoll(path_aniso_B2C, linearGradient), \
#                  computeaboveRoll(path_aniso_C2B_02, linearGradient) + computeaboveRoll(path_aniso_B2C_02, linearGradient), \
#                  facecolor = 'b',\
#                  alpha = 0.5)
#
#a2 = ax.fill_between(linearGradient, \
#                  computeaboveRoll(path_iso_C2B, linearGradient) + computeaboveRoll(path_iso_B2C, linearGradient), \
#                  computeaboveRoll(path_iso_C2B_02, linearGradient) + computeaboveRoll(path_iso_B2C_02, linearGradient), \
#                  facecolor = 'orange',\
#                  alpha = 0.5)
#pp = computeaboveRoll(path_aniso_B2C, linearGradient)
#pp[np.where(pp == 0)] = np.nan
#ax.plot(linearGradient,pp, 'b')
#pp = computeaboveRoll(path_aniso_B2C_02, linearGradient)
#pp[np.where(pp == 0)] = np.nan
#ax.plot(linearGradient, pp, 'b')
#pp = computeaboveRoll(path_iso_B2C, linearGradient)
#pp[np.where(pp == 0)] = np.nan
#ax.plot(linearGradient,pp, 'orange')
#pp = computeaboveRoll(path_iso_B2C_02, linearGradient)
#pp[np.where(pp == 0)] = np.nan
#ax.plot(linearGradient, pp, 'orange')
#ax.legend((a1,a2),('B2C - Anisotropic', 'B2C - Isotropic'))
#ax.set_xlabel('Absolute Roll Threshold [degrees]')
#ax.set_ylabel('Time above absolute roll threshold [m]')
#ax.set_xlim([0,16])

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
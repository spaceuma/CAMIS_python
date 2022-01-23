"""
MIT License
-----------

Copyright (c) 2021 University of Malaga
Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

Authors: J. Ricardo Sánchez Ibáñez, Carlos J. Pérez del Pulgar
Affiliation: Department of Systems Engineering and Automation
Space Robotics Lab (www.uma.es/space-robotics)
"""

# -*- coding: utf-8 -*-

import numpy as np

import data.cuadrigaData.cuadriga_reader as cr
from context import camis
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yaml
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import copy
import statistics
import scipy.signal
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy.spatial.transform import Rotation as R
from scipy import interpolate
import plotly.graph_objects as go
from matplotlib.colors import LightSource

from scipy.interpolate import interp1d

## =============================================================================
### IMPORT PLANNED PATHS
## =============================================================================
# Import 20201030execution.spydata using the variable explorer
#hiRes_elevationMap = np.loadtxt(\
#                        open("data/terrainData/UMATerrainCuesta_10cmDEM.csv",\
#                        "rb"), delimiter=" ", skiprows=0)
#    
#hiRes = 0.1
#offset = np.loadtxt(open("data/terrainData/UMATerrainCuesta_10cmOffset.csv",\
#                                 "rb"), delimiter=" ", skiprows=0)
#
#
#
#hexRes = 1.0/np.sqrt(6)
#
#localDEM = hiRes_elevationMap[260:510,350:600]
#offset = offset + [350*0.1, 260*0.1]
#occupancy_radius = 0.5
#tracking_error = 0.5
#env = camis.AnisotropicMap(hiRes_elevationMap[260:510,350:600], hiRes, hexRes,\
#                               offset, occupancy_radius, tracking_error)
#

#plt.style.use('default')
#plt.rcParams["font.family"] = "Constantia"
#plt.rcParams['mathtext.fontset'] = 'cm'
#plt.rcParams['mathtext.rm'] = 'serif'
#fig, ax1 = plt.subplots(figsize=(3.2, 3.9),nrows = 1, ncols = 1, constrained_layout=True)
#cc = ax1.scatter(env.xMap, env.yMap, c = 180/np.pi*np.arctan2(env.aspectY,env.aspectX),
#                 cmap="gist_rainbow",s=16.0, vmin = -180.0, vmax = 180.0, rasterized=True)
#cbar = fig.colorbar(cc, ax=ax1,shrink=0.6, location = 'top', ticks=[-180, -90, 0, 90, 180])
#cbar.set_label('Aspect Angle [deg]')
#ax1.set_xlim([0,24.5])
#ax1.set_ylim([0,24.5])
#ax1.set_xlabel('X [m]')
#ax1.set_ylabel('Y [m]')
#ax1.set_aspect('equal')
            
            
def buildPathInfo(posX, posY, time, Current, Heading):
    posX = posX - env.offset[0]
    posY = posY - env.offset[1]
    time = time - time[0]
    
    
    
    #posX, posY, time, roll, pitch, yaw, current, segment, heading, beta, slope
    pathInfo = np.concatenate(( posX, \
                                posY, \
                                time, time, time, time, \
                                Current, time, Heading, time, time)).reshape((11, len(posX)))
    
    
    for index, waypoint in enumerate(pathInfo[0]):
        # Segment length and heading angle
        if index == 0:
            A = [pathInfo[0][1] - pathInfo[0][0], pathInfo[1][1] - pathInfo[1][0]]
            pathInfo[7][index] = np.linalg.norm(A) * 0.5
#            pathInfo[8][index] = np.arctan2(A[1],A[0]) * 180.0/3.1416
        elif index == pathInfo[0].size - 1:
            A = [pathInfo[0][index] - pathInfo[0][index-1], pathInfo[1][index] - pathInfo[1][index-1]]
            pathInfo[7][index] = np.linalg.norm(A) * 0.5
#            pathInfo[8][index] = np.arctan2(A[1],A[0]) * 180.0/3.1416
        else:
            A1 = [pathInfo[0][index] - pathInfo[0][index-1], pathInfo[1][index] - pathInfo[1][index-1]]
            A2 = [pathInfo[0][index+1] - pathInfo[0][index], pathInfo[1][index+1] - pathInfo[1][index]]
            pathInfo[7][index] = np.linalg.norm(A1) * 0.5 + np.linalg.norm(A2) * 0.5
#            pathInfo[8][index] = np.arctan2(A1[1]+A2[1],A1[0]+A1[0]) * 180.0/3.1416

    
    B = env.getBeta(posX, posY, pathInfo[8])
    G = env.getSlope(posX, posY, pathInfo[8])
    betaX = np.cos(B)
    betaY = np.sin(B)
    B = np.arctan2(betaY,betaX)
    G = np.abs(G)*180/np.pi
    
    for i,_ in enumerate(G):
        pathInfo[4][i] = camis.ab2pitch(G[i],B[i])
        if B[i] >  np.pi/2 or B[i] < -np.pi/2:
            pathInfo[4][i] = - pathInfo[4][i]
        pathInfo[9][i] = B[i]
        pathInfo[10][i] = G[i]
        pathInfo[8][i] = pathInfo[8][i]*180/np.pi
    
    return pathInfo



### A2D results
posX_A2D, posY_A2D, Heading_A2D, Roll_A2D, Pitch_A2D, Yaw_A2D, \
Current_A2D, Speed_A2D, Distance_A2D, Segment_A2D, _, dTime_A2D, Time_A2D, \
Slope_A2D = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_06_43.txt',2)
 
posX_A2D_02, posY_A2D_02, Heading_A2D_02, Roll_A2D_02, Pitch_A2D_02, Yaw_A2D_02, Current_A2D_02, _, Distance_A2D, Segment_A2D, _, \
 dTime_A2D, Time_A2D_02, Slope_A2D_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_19_51.txt',1)

path_aniso_A2D = buildPathInfo(posX_A2D, posY_A2D, Time_A2D, Current_A2D,\
                               Heading_A2D)
path_aniso_A2D_02 = buildPathInfo(posX_A2D_02, posY_A2D_02, Time_A2D_02, \
                                  Current_A2D_02,Heading_A2D_02)


### D2A results
posX_D2A, posY_D2A, Heading_D2A, Roll_D2A, Pitch_D2A, Yaw_D2A, Current_D2A, _, Distance_D2A, Segment_D2A, _, \
 dTime_D2A, Time_D2A, Slope_D2A = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_10_23.txt',1)
posX_D2A_02, posY_D2A_02, Heading_D2A_02, Roll_D2A_02, Pitch_D2A_02, Yaw_D2A_02, Current_D2A_02, _, Distance_D2A, Segment_D2A, _, \
 dTime_D2A, Time_D2A_02, Slope_D2A_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_22_22.txt',1)
path_aniso_D2A = buildPathInfo(posX_D2A, posY_D2A, Time_D2A, Current_D2A, \
                               Heading_D2A)
path_aniso_D2A_02 = buildPathInfo(posX_D2A_02, posY_D2A_02, Time_D2A_02, \
                                  Current_D2A_02, Heading_D2A_02)


### C2B Results
posX_C2B, posY_C2B, Heading_C2B, Roll_C2B, Pitch_C2B, Yaw_C2B, \
Current_C2B, Speed, Distance_C2B, Segment_C2B, GPSspeed, \
 dTime_C2B, Time_C2B, Slope_C2B = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_31_06.txt',1)
posX_C2B_02, posY_C2B_02, Heading_C2B_02, Roll_C2B_02, Pitch_C2B_02, Yaw_C2B_02, \
Current_C2B_02, Speed, Distance_C2B, Segment_C2B, GPSspeed, \
 dTime_C2B, Time_C2B_02, Slope_C2B_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_41_07.txt',1)
path_aniso_C2B = buildPathInfo(posX_C2B, posY_C2B, Time_C2B, Current_C2B, \
                               Heading_C2B)
path_aniso_C2B_02 = buildPathInfo(posX_C2B_02, posY_C2B_02, Time_C2B_02, \
                                  Current_C2B_02, Heading_C2B_02)
path_aniso_C2B_02[2] = path_aniso_C2B_02[2] - 6


### B2C Results
posX_B2C, posY_B2C, Heading_B2C, Roll_B2C, Pitch_B2C, Yaw_B2C, \
Current_B2C, Speed, Distance_B2C, Segment_B2C, GPSspeed, \
 dTime_B2C, Time_B2C, Slope_B2C = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_33_17.txt',1)
posX_B2C_02, posY_B2C_02, Heading_B2C_02, Roll_B2C_02, Pitch_B2C_02, Yaw_B2C_02, \
Current_B2C_02, Speed, Distance_B2C, Segment_B2C, GPSspeed, \
 dTime_B2C, Time_B2C_02, Slope_B2C_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_43_28.txt',1)
path_aniso_B2C = buildPathInfo(posX_B2C, posY_B2C, Time_B2C, Current_B2C,\
                               Heading_B2C)
path_aniso_B2C_02 = buildPathInfo(posX_B2C_02, posY_B2C_02, Time_B2C_02, \
                                  Current_B2C_02, Heading_B2C_02)







fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
ax.plot(path_aniso_A2D[2], path_aniso_A2D[4])
#ax.plot(Time_A2D,[x*180/np.pi for x in Heading_A2D])



#### CUADRIGA COMPARATIVE ####
plannedDistA2D = np.cumsum(env_CUAD03_scene01[4].pathSegment)*2
plannedDistD2A = np.cumsum(env_CUAD03_scene01[3].pathSegment)*2
plannedDistB2C = np.cumsum(env_CUAD03_scene01[1].pathSegment)*2
plannedDistC2B = np.cumsum(env_CUAD03_scene01[6].pathSegment)*2
path_aniso_A2D[2] = path_aniso_A2D[2] - path_aniso_A2D[2][-1] + \
                    plannedDistA2D[-1] - 2.0 #2.0 = lookahead?
path_aniso_A2D_02[2] = path_aniso_A2D_02[2] - path_aniso_A2D_02[2][-1] + \
                       plannedDistA2D[-1] - 1.0 #2.0 = lookahead?
#path_iso_A2D[2] = path_iso_A2D[2] - path_iso_A2D[2][-1] + \
#                    plannedDistA2D[-1] - 2.0 #2.0 = lookahead?
#path_iso_A2D_02[2] = path_iso_A2D_02[2] - path_iso_A2D_02[2][-1] + \
#                       plannedDistA2D[-1] - 2.0 #2.0 = lookahead?
path_aniso_D2A[2] = path_aniso_D2A[2] - path_aniso_D2A[2][-1] + \
                    plannedDistD2A[-1] - 2.0 #2.0 = lookahead?
path_aniso_D2A_02[2] = path_aniso_D2A_02[2] - path_aniso_D2A_02[2][-1] + \
                       plannedDistD2A[-1] - 2.0 #2.0 = lookahead?
path_aniso_B2C[2] = path_aniso_B2C[2] - path_aniso_B2C[2][-1] + \
                    plannedDistB2C[-1] - 4.0#2.0 = lookahead?
path_aniso_B2C_02[2] = path_aniso_B2C_02[2] - path_aniso_B2C_02[2][-1] + \
                       plannedDistB2C[-1] - 2.0 #2.0 = lookahead?
path_aniso_C2B[2] = path_aniso_C2B[2] - path_aniso_C2B[2][-1] + \
                    plannedDistC2B[-1] - 2.0 #2.0 = lookahead?
path_aniso_C2B_02[2] = path_aniso_C2B_02[2] - path_aniso_C2B_02[2][-1] + \
                       plannedDistC2B[-1] - 2.0 #2.0 = lookahead?


plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, rowaxes = plt.subplots(figsize=(5, 6), nrows = 4, ncols = 1)
plt.subplots_adjust(left = 0.085, right = 0.95, bottom = 0.075, \
                    top = 0.98, wspace = 0.1, hspace = 0.1)
fig2, rowaxes2 = plt.subplots(figsize=(5, 6), nrows = 4, ncols = 1)
plt.subplots_adjust(left = 0.085, right = 0.95, bottom = 0.075, \
                    top = 0.98, wspace = 0.1, hspace = 0.1)

deg2rad = np.pi / 180.0

timeArray = []
costArray = []
pitchArray = []
rollArray = []
pitchArrayFiltered = []
rollArrayFiltered = []

for i,path in enumerate((path_aniso_A2D, path_aniso_A2D_02,\
                        path_aniso_D2A, path_aniso_D2A_02,\
                        path_aniso_B2C, path_aniso_B2C_02, \
                        path_aniso_C2B, path_aniso_C2B_02)):
 
    p1 = rowaxes[(int)(i/2)].plot(path[2], path[4], 'c', alpha = 0.2)
    rollarray = np.arccos(np.cos(path[10]*deg2rad)/np.cos(path[4]*deg2rad))*180.0 / np.pi
    p2 = rowaxes[(int)(i/2)].plot(path[2], rollarray, 'r', alpha = 0.2)
    
    if i%2 == 0:
        timeArray = path[2]
        pitchArray = path[4]
        rollArray = rollarray
    else:
        timeArray = np.concatenate((timeArray, path[2]))
        pitchArray = np.concatenate((pitchArray, path[4]))
        rollArray = np.concatenate((rollArray, rollarray))
        sortedPitchArray = [x for _,x in sorted(zip(timeArray,pitchArray))]
        sortedRollArray = [x for _,x in sorted(zip(timeArray,rollArray))]
        
        sortedPitchArray = [0 if x != x else x for x in sortedPitchArray]
        sortedRollArray = [0 if x != x else x for x in sortedRollArray]
        pitchArrayFiltered = savgol_filter(sortedPitchArray, 51, 3)
        rollArrayFiltered = savgol_filter(sortedRollArray, 51, 3)
        p1 = rowaxes[(int)(i/2)].plot(sorted(timeArray), pitchArrayFiltered, 'c')
        p2 = rowaxes[(int)(i/2)].plot(sorted(timeArray), rollArrayFiltered, 'r')
    
    rowaxes2[(int)(i/2)].plot(path[2], \
        path[6]*np.append(0,np.diff(path[2]))/path[7], 'orange', alpha = 0.2)
    
    if i%2 == 0:
        timeArray = path[2]
        costArray = path[6]*np.append(0,np.diff(path[2]))/path[7]
    else:
#        timeArray = np.concatenate((timeArray, path[2]))
        costArray = np.concatenate((costArray, path[6]*np.append(0,np.diff(path[2]))/path[7]))
        sortedCostArray = [x for _,x in sorted(zip(timeArray,costArray))]
        costArrayFiltered = savgol_filter(sortedCostArray, 51, 3)
        p7, = rowaxes2[(int)(i/2)].plot(sorted(timeArray), costArrayFiltered, 'orange')
        
    if i==1:
        pitchArrayFilteredAD = pitchArrayFiltered
        rollArrayFilteredAD = rollArrayFiltered
        costArrayFilteredAD = costArrayFiltered
        timeAD = sorted(timeArray)
    if i==3:
        pitchArrayFilteredDA = pitchArrayFiltered
        rollArrayFilteredDA = rollArrayFiltered
        costArrayFilteredDA = costArrayFiltered
        timeDA = sorted(timeArray)
    if i==5:
        pitchArrayFilteredBC = pitchArrayFiltered
        rollArrayFilteredBC = rollArrayFiltered
        costArrayFilteredBC = costArrayFiltered
        timeBC = sorted(timeArray)
    if i==7:
        pitchArrayFilteredCB = pitchArrayFiltered
        rollArrayFilteredCB = rollArrayFiltered
        costArrayFilteredCB = costArrayFiltered
        timeCB = sorted(timeArray)
        
        
    #
    

rowaxes[3].set_xlabel('Elapsed Time [s]')
rowaxes2[3].set_xlabel('Elapsed Time [s]')
fig.text(0.01, 0.5, 'Orientation Angle [deg]', va='center', rotation='vertical')
fig2.text(0.01, 0.5, 'Energy per meter [As/m]', va='center', rotation='vertical')

p3 = rowaxes[0].plot(plannedDistA2D, \
        180.0/np.pi*np.asarray(env_CUAD03_scene01[4].pathRoll), 'orange')
rowaxes[0].text(49.0, 6.0, '$\mathbf{x_a \Rightarrow x_d}$', ha='right', fontsize=12)
rowaxes[1].plot(plannedDistD2A, \
        180.0/np.pi*np.asarray(env_CUAD03_scene01[3].pathRoll), 'orange')
rowaxes[1].text(49.0, 16.0, '$\mathbf{x_d \Rightarrow x_a}$', ha='right', fontsize=12)
rowaxes[2].plot(plannedDistB2C, \
        180.0/np.pi*np.asarray(env_CUAD03_scene01[1].pathRoll), 'orange')
rowaxes[2].text(49.0, 15.0, '$\mathbf{x_b \Rightarrow x_c}$', ha='right', fontsize=12)
rowaxes[3].plot(plannedDistC2B, \
        180.0/np.pi*np.asarray(env_CUAD03_scene01[6].pathRoll), 'orange')
rowaxes[3].text(49.0, 5.0, '$\mathbf{x_c \Rightarrow x_b}$', ha='right', fontsize=12)
p4 = rowaxes[0].plot(plannedDistA2D, \
        180.0/np.pi*np.asarray(env_CUAD03_scene01[4].pathPitch), 'b')
rowaxes[1].plot(plannedDistD2A, \
        180.0/np.pi*np.asarray(env_CUAD03_scene01[3].pathPitch), 'b')
rowaxes[2].plot(plannedDistB2C, \
        180.0/np.pi*np.asarray(env_CUAD03_scene01[1].pathPitch), 'b')
rowaxes[3].plot(plannedDistC2B, \
        180.0/np.pi*np.asarray(env_CUAD03_scene01[6].pathPitch), 'b')

p6, = rowaxes2[0].plot(plannedDistA2D, \
        np.asarray(env_CUAD03_scene01[4].pathCost)*0.5, 'b')
rowaxes2[0].text(1.0, 140.0, '$\mathbf{x_a \Rightarrow x_d}$', ha='left', fontsize=12)
rowaxes2[1].plot(plannedDistD2A, \
        np.asarray(env_CUAD03_scene01[3].pathCost)*0.5, 'b')
rowaxes2[1].text(1.0, 140.0, '$\mathbf{x_d \Rightarrow x_a}$', ha='left', fontsize=12)
rowaxes2[2].plot(plannedDistB2C, \
        np.asarray(env_CUAD03_scene01[1].pathCost)*0.5, 'b')
rowaxes2[2].text(1.0, 90.0, '$\mathbf{x_b \Rightarrow x_c}$', ha='left', fontsize=12)
rowaxes2[3].plot(plannedDistC2B, \
        np.asarray(env_CUAD03_scene01[6].pathCost)*0.5, 'b')
rowaxes2[3].text(1.0, 70.0, '$\mathbf{x_c \Rightarrow x_b}$', ha='left', fontsize=12)
for ax in rowaxes:
    ax.set_xlim([0,50])
    ax.tick_params(axis="x",direction="in", pad=-3)
for ax in rowaxes2:
    ax.set_xlim([0,50])
    ax.tick_params(axis="x",direction="in", pad=-3)
p5, = ax.plot([0], marker='None',
           linestyle='None', label='Wheel Model')
rowaxes[0].legend((p5,p2[0],p1[0],p5,p3[0],p4[0]),\
       ('Experimental','Roll','Pitch',\
        'Planned', 'Roll','Pitch'), ncol = 2)
rowaxes2[0].legend((p7,p6), ('Experimental','Planned'))





def getStats(plan_time, plan_data, exp_time, exp_data):
    f = interp1d(plan_time, plan_data, kind='cubic')
    xdata1 = [item for item in exp_time if item >= plan_time[0]] 
    ydata1 = exp_data[len(exp_time)-len(xdata1)-1:-1]
    ydata2 = f(xdata1)
    RMSE = 0
    MAE = 0
    for i,p in enumerate(xdata1):
        RMSE = RMSE + (ydata1[i] - ydata2[i])**2
        MAE = MAE + np.abs(ydata1[i] - ydata2[i])
    RMSE_raw = RMSE
    MAE_raw = MAE
    N = len(ydata1)
    RMSE = np.sqrt(RMSE / len(ydata1))
    MAE = MAE / len(ydata1)
    return RMSE, MAE, RMSE_raw, MAE_raw, N

def getResults(plan_time, enviro, exp_time, exp_roll, exp_pitch, exp_cost):
    RMSE_roll, MAE_roll, RMSE_raw_roll, MAE_raw_roll, N_roll = getStats(plan_time,
                                   180.0/np.pi*np.asarray(enviro.pathRoll),
                                   exp_time,
                                   exp_roll)
    RMSE_pitch, MAE_pitch, RMSE_raw_pitch, MAE_raw_pitch, N_pitch = getStats(plan_time,
                                     180.0/np.pi*np.asarray(enviro.pathPitch),
                                     exp_time,
                                     exp_pitch)
    RMSE_cost, MAE_cost, RMSE_raw_cost, MAE_raw_cost, N_cost = getStats(plan_time,
                                     np.asarray(enviro.pathCost)*0.5,
                                     exp_time,
                                     exp_cost)
    print('Roll RMSE = ')
    print(RMSE_roll)
    print('Roll MAE = ')
    print(MAE_roll)
    print('Pitch RMSE = ')
    print(RMSE_pitch)
    print('Pitch MAE = ')
    print(MAE_pitch)
    print('Cost RMSE = ')
    print(RMSE_cost)
    print('Cost MAE = ')
    print(MAE_cost)
    
    return RMSE_raw_roll, MAE_raw_roll, N_roll, RMSE_raw_pitch, MAE_raw_pitch, \
           N_pitch,RMSE_raw_cost, MAE_raw_cost, N_cost

RMSE_raw_rollAD, MAE_raw_rollAD, N_rollAD, RMSE_raw_pitchAD, MAE_raw_pitchAD, \
N_pitchAD, RMSE_raw_costAD, MAE_raw_costAD, N_costAD = \
getResults(plannedDistA2D, env_CUAD03_scene01[4], timeAD, rollArrayFilteredAD,
           pitchArrayFilteredAD, costArrayFilteredAD)

RMSE_raw_rollDA, MAE_raw_rollDA, N_rollDA, RMSE_raw_pitchDA, MAE_raw_pitchDA, \
           N_pitchDA, RMSE_raw_costDA, MAE_raw_costDA, N_costDA = \
getResults(plannedDistD2A, env_CUAD03_scene01[3], timeDA, rollArrayFilteredDA,
           pitchArrayFilteredDA, costArrayFilteredDA)

RMSE_raw_rollBC, MAE_raw_rollBC, N_rollBC, RMSE_raw_pitchBC, MAE_raw_pitchBC, \
           N_pitchBC, RMSE_raw_costBC, MAE_raw_costBC, N_costBC = \
getResults(plannedDistB2C, env_CUAD03_scene01[1], timeBC, rollArrayFilteredBC,
           pitchArrayFilteredBC, costArrayFilteredBC)

RMSE_raw_rollCB, MAE_raw_rollCB, N_rollCB, RMSE_raw_pitchCB, MAE_raw_pitchCB, \
           N_pitchCB, RMSE_raw_costCB, MAE_raw_costCB, N_costCB = \
getResults(plannedDistC2B, env_CUAD03_scene01[6], timeCB, rollArrayFilteredCB,
           pitchArrayFilteredCB, costArrayFilteredCB)


RMSE_roll_totalAD = RMSE_raw_rollAD + RMSE_raw_rollDA + RMSE_raw_rollBC + RMSE_raw_rollCB
RMSE_roll_totalAD = np.sqrt(RMSE_roll_totalAD / (N_rollAD + N_rollDA + N_rollBC + N_rollCB))
MAE_roll_totalAD = MAE_raw_rollAD + MAE_raw_rollDA + MAE_raw_rollBC + MAE_raw_rollCB
MAE_roll_totalAD = MAE_roll_totalAD / (N_rollAD + N_rollDA + N_rollBC + N_rollCB)

RMSE_pitch_totalAD = RMSE_raw_pitchAD + RMSE_raw_pitchDA + RMSE_raw_pitchBC + RMSE_raw_pitchCB
RMSE_pitch_totalAD = np.sqrt(RMSE_pitch_totalAD / (N_pitchAD + N_pitchDA + N_pitchBC + N_pitchCB))
MAE_pitch_totalAD = MAE_raw_pitchAD + MAE_raw_pitchDA + MAE_raw_pitchBC + MAE_raw_pitchCB
MAE_pitch_totalAD = MAE_pitch_totalAD / (N_pitchAD + N_pitchDA + N_pitchBC + N_pitchCB)

RMSE_cost_totalAD = RMSE_raw_costAD + RMSE_raw_costDA + RMSE_raw_costBC + RMSE_raw_costCB
RMSE_cost_totalAD = np.sqrt(RMSE_cost_totalAD / (N_costAD + N_costDA + N_costBC + N_costCB))
MAE_cost_totalAD = MAE_raw_costAD + MAE_raw_costDA + MAE_raw_costBC + MAE_raw_costCB
MAE_cost_totalAD = MAE_cost_totalAD / (N_costAD + N_costDA + N_costBC + N_costCB)

print('Roll RMSE = ')
print(RMSE_roll_totalAD)
print('Roll MAE = ')
print(MAE_roll_totalAD)
print('Pitch RMSE = ')
print(RMSE_pitch_totalAD)
print('Pitch MAE = ')
print(MAE_pitch_totalAD)
print('Cost RMSE = ')
print(RMSE_cost_totalAD)
print('Cost MAE = ')
print(MAE_cost_totalAD)



RMSE_AD = 0
MAE_AD = 0
for i,p in enumerate(pitchAD_exp):
    RMSE_AD = RMSE_AD + (pitchAD_exp[i] - pitchAD_planned[i])**2
    MAE_AD = MAE_AD + np.abs(pitchAD_exp[i] - pitchAD_planned[i])
RMSE_AD = np.sqrt(RMSE_AD / len(timeAD_exp))
MAE_AD = MAE_AD / len(timeAD_exp)


fig_comp, axes_comp = plt.subplots()

axes_comp.plot(timeAD_exp,pitchAD_exp)
axes_comp.plot(timeAD_exp, pitchAD_planned)

### ISOTROPIC CASE

### A2D results
posX_isoA2D, posY_isoA2D, Heading_isoA2D, Roll_isoA2D, Pitch_isoA2D, Yaw_A2D, \
Current_isoA2D, Speed_A2D, Distance_A2D, Segment_A2D, GPSspeed_A2D, \
 dTime_isoA2D, Time_isoA2D, Slope_isoA2D = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_13_37.txt',1)
posX_isoA2D_02, posY_isoA2D_02, Heading_isoA2D_02, Roll_isoA2D_02, Pitch_isoA2D_02, Yaw_A2D_02, \
Current_isoA2D_02, Speed_isoA2D_02, Distance_isoA2D_02, Segment_isoA2D_02, GPSspeed_isoA2D_02, \
 dTime_isoA2D, Time_isoA2D_02, Slope_isoA2D_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_24_59.txt',1)
path_iso_A2D = buildPathInfo(posX_isoA2D, posY_isoA2D, Time_isoA2D, \
                             Current_isoA2D, Heading_isoA2D)
path_iso_A2D_02 = buildPathInfo(posX_isoA2D_02, posY_isoA2D_02, Time_isoA2D_02, \
                                Current_isoA2D_02, Heading_isoA2D_02)


### D2A results
posX_isoD2A, posY_isoD2A, Heading_isoD2A, Roll_isoD2A, Pitch_isoD2A, Yaw_isoD2A, \
Current_isoD2A, Speed_D2A, Distance_D2A, Segment_D2A, GPSspeed, \
 dTime_D2A, Time_isoD2A, Slope_isoD2A = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_27_26.txt',1)
posX_isoD2A_02, posY_isoD2A_02, Heading_isoD2A_02, Roll_isoD2A_02, Pitch_isoD2A_02, Yaw_isoD2A_02, \
Current_isoD2A_02, Speed_D2A, Distance_D2A, Segment_D2A, GPSspeed, \
 dTime_D2A, Time_isoD2A_02, Slope_isoD2A_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_12_06_16.txt',1)
path_iso_D2A_02 = buildPathInfo(posX_isoD2A_02, posY_isoD2A_02, Time_isoD2A_02, \
                                Current_isoD2A_02, Heading_isoD2A_02)
path_iso_D2A = buildPathInfo(posX_isoD2A, posY_isoD2A, Time_isoD2A, \
                             Current_isoD2A, Heading_isoD2A)


### C2B results
posX_isoC2B, posY_isoC2B, Heading_isoC2B, Roll_isoC2B, Pitch_isoC2B, Yaw_isoC2B, \
Current_isoC2B, Speed, Distance_C2B, Segment_C2B, GPSspeed, \
 dTime_C2B, Time_isoC2B, Slope_isoC2B = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_35_30.txt',1)
posX_isoC2B_02, posY_isoC2B_02, Heading_isoC2B_02, Roll_isoC2B_02, Pitch_isoC2B_02, Yaw_isoC2B_02, \
Current_isoC2B_02, Speed, Distance_C2B, Segment_C2B, GPSspeed, \
 dTime_C2B, Time_isoC2B_02, Slope_isoC2B_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_45_47.txt',1)
path_iso_C2B = buildPathInfo(posX_isoC2B, posY_isoC2B, Time_isoC2B, \
                             Current_isoC2B, Heading_isoC2B)
path_iso_C2B_02 = buildPathInfo(posX_isoC2B_02, posY_isoC2B_02, Time_isoC2B_02, \
                                Current_isoC2B_02, Heading_isoC2B_02)


### B2C results
posX_isoB2C, posY_isoB2C, Heading_isoB2C, Roll_isoB2C, Pitch_isoB2C, Yaw_isoB2C, \
Current_isoB2C, Speed, Distance_B2C, Segment_B2C, GPSspeed, \
 dTime_isoB2C, Time_isoB2C, Slope_isoB2C = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_38_48.txt',1)
posX_isoB2C_02, posY_isoB2C_02, Heading_isoB2C_02, Roll_isoB2C_02, Pitch_isoB2C_02, Yaw_isoB2C_02, \
Current_isoB2C_02, Speed, Distance_B2C, Segment_B2C, GPSspeed, \
 dTime_B2C, Time_isoB2C_02, Slope_isoB2C_02 = cr.readCuadrigaData('experimental_results/2020_11_03/2020_11_03_11_48_42.txt',1)
path_iso_B2C = buildPathInfo(posX_isoB2C, posY_isoB2C, Time_isoB2C, \
                             Current_isoB2C, Heading_isoB2C)
path_iso_B2C_02 = buildPathInfo(posX_isoB2C_02, posY_isoB2C_02, Time_isoB2C_02, \
                                Current_isoB2C_02, Heading_isoB2C_02) 





##### CUADRIGA COMPARATIVE ####
#plannedDistA2D = np.cumsum(env_isoCUAD03_scene01[4].pathSegment)*2
#plannedDistD2A = np.cumsum(env_isoCUAD03_scene01[3].pathSegment)*2
#plannedDistB2C = np.cumsum(env_isoCUAD03_scene01[1].pathSegment)*2
#plannedDistC2B = np.cumsum(env_isoCUAD03_scene01[6].pathSegment)*2
##path_aniso_A2D[2] = path_aniso_A2D[2] - path_aniso_A2D[2][-1] + \
##                    plannedDistA2D[-1] - 2.0 #2.0 = lookahead?
##path_aniso_A2D_02[2] = path_aniso_A2D_02[2] - path_aniso_A2D_02[2][-1] + \
##                       plannedDistA2D[-1] - 1.0 #2.0 = lookahead?
#path_iso_A2D[2] = path_iso_A2D[2] - path_iso_A2D[2][-1] + \
#                    plannedDistA2D[-1] - 4.0 #2.0 = lookahead?
#path_iso_A2D_02[2] = path_iso_A2D_02[2] - path_iso_A2D_02[2][-1] + \
#                       plannedDistA2D[-1] - 4.0 #2.0 = lookahead?
#path_iso_D2A[2] = path_iso_D2A[2] - path_iso_D2A[2][-1] + \
#                    plannedDistD2A[-1] - 0.0 #2.0 = lookahead?
#path_iso_D2A_02[2] = path_iso_D2A_02[2] - path_iso_D2A_02[2][-1] + \
#                       plannedDistD2A[-1] - 0.0 #2.0 = lookahead?
#path_iso_C2B[2] = path_iso_C2B[2] - path_iso_C2B[2][-1] + \
#                    plannedDistC2B[-1] + 0.0 #2.0 = lookahead?
#path_iso_C2B_02[2] = path_iso_C2B_02[2] - path_iso_C2B_02[2][-1] + \
#                       plannedDistC2B[-1] + 0.0 #2.0 = lookahead?                     
#path_iso_B2C[2] = path_iso_B2C[2] - path_iso_B2C[2][-1] + \
#                    plannedDistB2C[-1] - 0.0 #2.0 = lookahead?
#path_iso_B2C_02[2] = path_iso_B2C_02[2] - path_iso_B2C_02[2][-1] + \
#                       plannedDistB2C[-1] - 0.0 #2.0 = lookahead?                     
#                                              
#                       
##path_aniso_D2A[2] = path_aniso_D2A[2] - path_aniso_D2A[2][-1] + \
##                    plannedDistD2A[-1] - 2.0 #2.0 = lookahead?
##path_aniso_D2A_02[2] = path_aniso_D2A_02[2] - path_aniso_D2A_02[2][-1] + \
##                       plannedDistD2A[-1] - 2.0 #2.0 = lookahead?
##path_aniso_B2C[2] = path_aniso_B2C[2] - path_aniso_B2C[2][-1] + \
##                    plannedDistB2C[-1] - 4.0#2.0 = lookahead?
##path_aniso_B2C_02[2] = path_aniso_B2C_02[2] - path_aniso_B2C_02[2][-1] + \
##                       plannedDistB2C[-1] - 2.0 #2.0 = lookahead?
##path_aniso_C2B[2] = path_aniso_C2B[2] - path_aniso_C2B[2][-1] + \
##                    plannedDistC2B[-1] - 2.0 #2.0 = lookahead?
##path_aniso_C2B_02[2] = path_aniso_C2B_02[2] - path_aniso_C2B_02[2][-1] + \
##                       plannedDistC2B[-1] - 2.0 #2.0 = lookahead?
#
#
#plt.style.use('seaborn-darkgrid')
#plt.rcParams["font.family"] = "Constantia"
#plt.rcParams['mathtext.fontset'] = 'cm'
#plt.rcParams['mathtext.rm'] = 'serif'
#fig, rowaxes = plt.subplots(figsize=(5, 6), nrows = 4, ncols = 1)
#plt.subplots_adjust(left = 0.085, right = 0.95, bottom = 0.075, \
#                    top = 0.98, wspace = 0.1, hspace = 0.1)
#fig2, rowaxes2 = plt.subplots(figsize=(5, 6), nrows = 4, ncols = 1)
#plt.subplots_adjust(left = 0.085, right = 0.95, bottom = 0.075, \
#                    top = 0.98, wspace = 0.1, hspace = 0.1)
#
#deg2rad = np.pi / 180.0
#
#timeArray = []
#costArray = []
#pitchArray = []
#rollArray = []
#pitchArrayFiltered = []
#rollArrayFiltered = []
#
#for i,path in enumerate((path_iso_A2D, path_iso_A2D_02,\
#                        path_iso_D2A, path_iso_D2A_02,\
#                        path_iso_B2C, path_iso_B2C_02, \
#                        path_iso_C2B, path_iso_C2B_02)):
# 
#    p1 = rowaxes[(int)(i/2)].plot(path[2], path[4], 'c', alpha = 0.2)
#    rollarray = np.arccos(np.cos(path[10]*deg2rad)/np.cos(path[4]*deg2rad))*180.0 / np.pi
#    p2 = rowaxes[(int)(i/2)].plot(path[2], rollarray, 'r', alpha = 0.2)
#    
#    if i%2 == 0:
#        timeArray = path[2]
#        pitchArray = path[4]
#        rollArray = rollarray
#    else:
#        timeArray = np.concatenate((timeArray, path[2]))
#        pitchArray = np.concatenate((pitchArray, path[4]))
#        rollArray = np.concatenate((rollArray, rollarray))
#        sortedPitchArray = [x for _,x in sorted(zip(timeArray,pitchArray))]
#        sortedRollArray = [x for _,x in sorted(zip(timeArray,rollArray))]        
#        sortedPitchArray = [0 if x != x else x for x in sortedPitchArray]
#        sortedRollArray = [0 if x != x else x for x in sortedRollArray]
#        pitchArrayFiltered = savgol_filter(sortedPitchArray, 51, 3)
#        rollArrayFiltered = savgol_filter(sortedRollArray, 51, 3)
#        p1 = rowaxes[(int)(i/2)].plot(sorted(timeArray), pitchArrayFiltered, 'c')
#        p2 = rowaxes[(int)(i/2)].plot(sorted(timeArray), rollArrayFiltered, 'r')
#    
#    
#    
#    
#    rowaxes2[(int)(i/2)].plot(path[2], \
#        path[6]*np.append(0,np.diff(path[2]))/path[7], 'orange', alpha = 0.2)
#    
#    if i%2 == 0:
#        timeArray = path[2]
#        costArray = path[6]*np.append(0,np.diff(path[2]))/path[7]
#    else:
##        timeArray = np.concatenate((timeArray, path[2]))
#        costArray = np.concatenate((costArray, path[6]*np.append(0,np.diff(path[2]))/path[7]))
#        sortedCostArray = [x for _,x in sorted(zip(timeArray,costArray))]
#        costArrayFiltered = savgol_filter(sortedCostArray, 51, 3)
#        p7, = rowaxes2[(int)(i/2)].plot(sorted(timeArray), costArrayFiltered, 'orange')
#    #
#    
#
#rowaxes[3].set_xlabel('Elapsed Time [s]')
#rowaxes2[3].set_xlabel('Elapsed Time [s]')
#fig.text(0.01, 0.5, 'Orientation Angle [deg]', va='center', rotation='vertical')
#fig2.text(0.01, 0.5, 'Energy per meter [As/m]', va='center', rotation='vertical')
#
#p3 = rowaxes[0].plot(plannedDistA2D, \
#        180.0/np.pi*np.asarray(env_isoCUAD03_scene01[4].pathRoll), 'orange')
#rowaxes[0].text(49.0, 6.0, '$\mathbf{x_a \Rightarrow x_d}$', ha='right', fontsize=12)
#rowaxes[1].plot(plannedDistD2A, \
#        180.0/np.pi*np.asarray(env_isoCUAD03_scene01[3].pathRoll), 'orange')
#rowaxes[1].text(49.0, 16.0, '$\mathbf{x_d \Rightarrow x_a}$', ha='right', fontsize=12)
#rowaxes[2].plot(plannedDistB2C, \
#        180.0/np.pi*np.asarray(env_isoCUAD03_scene01[1].pathRoll), 'orange')
#rowaxes[2].text(49.0, 15.0, '$\mathbf{x_b \Rightarrow x_c}$', ha='right', fontsize=12)
#rowaxes[3].plot(plannedDistC2B, \
#        180.0/np.pi*np.asarray(env_isoCUAD03_scene01[6].pathRoll), 'orange')
#rowaxes[3].text(49.0, 5.0, '$\mathbf{x_c \Rightarrow x_b}$', ha='right', fontsize=12)
#p4 = rowaxes[0].plot(plannedDistA2D, \
#        180.0/np.pi*np.asarray(env_isoCUAD03_scene01[4].pathPitch), 'b')
#rowaxes[1].plot(plannedDistD2A, \
#        180.0/np.pi*np.asarray(env_isoCUAD03_scene01[3].pathPitch), 'b')
#rowaxes[2].plot(plannedDistB2C, \
#        180.0/np.pi*np.asarray(env_isoCUAD03_scene01[1].pathPitch), 'b')
#rowaxes[3].plot(plannedDistC2B, \
#        180.0/np.pi*np.asarray(env_isoCUAD03_scene01[6].pathPitch), 'b')
#
#p6, = rowaxes2[0].plot(plannedDistA2D, \
#        np.asarray(env_isoCUAD03_scene01[4].pathCost)*0.5, 'b')
#rowaxes2[0].text(1.0, 140.0, '$\mathbf{x_a \Rightarrow x_d}$', ha='left', fontsize=12)
#rowaxes2[1].plot(plannedDistD2A, \
#        np.asarray(env_isoCUAD03_scene01[3].pathCost)*0.5, 'b')
#rowaxes2[1].text(1.0, 140.0, '$\mathbf{x_d \Rightarrow x_a}$', ha='left', fontsize=12)
#rowaxes2[2].plot(plannedDistB2C, \
#        np.asarray(env_isoCUAD03_scene01[1].pathCost)*0.5, 'b')
#rowaxes2[2].text(1.0, 90.0, '$\mathbf{x_b \Rightarrow x_c}$', ha='left', fontsize=12)
#rowaxes2[3].plot(plannedDistC2B, \
#        np.asarray(env_isoCUAD03_scene01[6].pathCost)*0.5, 'b')
#rowaxes2[3].text(1.0, 70.0, '$\mathbf{x_c \Rightarrow x_b}$', ha='left', fontsize=12)
#for ax in rowaxes:
#    ax.set_xlim([0,50])
#    ax.tick_params(axis="x",direction="in", pad=-3)
#for ax in rowaxes2:
#    ax.set_xlim([0,50])
#    ax.tick_params(axis="x",direction="in", pad=-3)
#p5, = ax.plot([0], marker='None',
#           linestyle='None', label='Wheel Model')
#rowaxes[0].legend((p5,p2[0],p1[0],p5,p3[0],p4[0]),\
#       ('Experimental','Roll','Pitch',\
#        'Planned', 'Roll','Pitch'), ncol = 2)
#rowaxes2[0].legend((p7,p6), ('Experimental','Planned'))


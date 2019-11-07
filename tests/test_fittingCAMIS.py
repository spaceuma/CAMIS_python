# -*- coding: utf-8 -*-
"""
DEMO script: building CAMIS based on Cuadriga experimental data
@author: J.Ricardo Sanchez Ibanez (ricardosan@uma.es), Github: @Ryk-San

This script is a reference demo to understand how CAMIS is built using data
obtained by a mobile robot.

"""

import numpy as np

import data.cuadrigaData.cuadriga_reader as cr
from context import camis
import matplotlib.pyplot as plt
import matplotlib.cm as cm


#fittingMode = 'DEM'
fittingMode = 'IMU'


if fittingMode == 'DEM':
    
    hiRes_elevationMap = np.loadtxt(open("data/terrainData/UMARescueArea_10cmDEM.csv",\
                                         "rb"), delimiter=" ", skiprows=0)
    hiRes_posX = np.loadtxt(open("data/terrainData/UMARescueArea_10cmPosX.csv",\
                                         "rb"), delimiter=" ", skiprows=0)
    hiRes_posY = np.loadtxt(open("data/terrainData/UMARescueArea_10cmPosY.csv",\
                                         "rb"), delimiter=" ", skiprows=0)
    
    hiRes = hiRes_posX[0,1] - hiRes_posX[0,0]
    env = camis.PDEM(hiRes_elevationMap, hiRes, .5, (hiRes_posX[0,0],hiRes_posY[0,0]))
    env.smoothMap(1.0) #Tracking error is not considered here



"01 - The CSV files collected during preliminar experiments are loaded"
csvFiles = [
            'data/cuadrigaData/20190624/2019_06_24_11_17_55.txt',
            'data/cuadrigaData/20190624/2019_06_24_11_22_21.txt',
            'data/cuadrigaData/20190624/2019_06_24_11_32_17.txt',
            'data/cuadrigaData/20190624/2019_06_24_11_36_06.txt',
            'data/cuadrigaData/20190624/2019_06_24_11_45_59.txt',
            'data/cuadrigaData/20190624/2019_06_24_11_49_29.txt',
            'data/cuadrigaData/20190624/2019_06_24_11_52_47.txt',
            'data/cuadrigaData/20190624/2019_06_24_11_57_35.txt'
            ]

#cr.showData('cuadrigaData/20190624/2019_06_24_11_17_55.txt')

"02 - CAMIS input data is extracted"
betaX = []
betaY = []
beta = []
cost = []
gradient = []
for file in csvFiles:
    " Cost C according to robot orientation Roll-Pitch-Yaw "
    posX, posY, heading, Roll, Pitch, Yaw, Current, Speed = cr.readCuadrigaData(file)
    " Conversion from Roll-Pitch-Yaw to Gradient-Beta_angle "
    if fittingMode == 'IMU':
        G, Bx, By = camis.rpy2ab(Roll,Pitch,Yaw)
        betaX = betaX + Bx
        betaY = betaY + By
        B = np.arctan2(By,Bx)
        beta = beta + B.tolist()
        currentCost = np.divide(Current,Speed)
        cost = cost + currentCost.tolist()
        gradient = gradient + G
    if fittingMode == 'DEM':
        G, B = env.getBeta(posX - hiRes_posX[0,0], posY - hiRes_posY[0,0], heading)
        betaX = np.cos(B)
        betaY = np.sin(B)
        B = np.arctan2(betaY,betaX)
        beta = beta + B.tolist()
        currentCost = np.divide(Current,Speed)
        cost = cost + currentCost.tolist()
        G = np.abs(G)*180/np.pi
        gradient = gradient + G.tolist()
    
    

if fittingMode == 'DEM':
    maxCost = max(cost)
    gradient_threshold = 20.0
if fittingMode == 'IMU':
    maxCost = 50.0
    gradient_threshold = max(gradient)
    

beta.append(np.pi)
cost.append(maxCost)
gradient.append(gradient_threshold)

beta.append(0)
cost.append(maxCost)
gradient.append(gradient_threshold)

beta.append(np.pi/2)
cost.append(maxCost)
gradient.append(gradient_threshold)

beta.append(-np.pi)
cost.append(maxCost)
gradient.append(gradient_threshold)



"03 - The functions Cd, Ca, Cl1 and Cl2 are obtained through a fitting process"
beta = [0 if np.isnan(x) else x for x in beta]
CdRoots, CaRoots, Cl1Roots, Cl2Roots = \
camis.computeDirCosts(gradient,beta,cost)
r1 = camis.CamisModel.fromRoots(CdRoots,CaRoots,Cl1Roots,Cl2Roots, gradient_threshold)

isoRoots = (np.asarray(CaRoots) + np.asarray(CdRoots) + np.asarray(Cl1Roots) + np.asarray(Cl2Roots))/4

rIsoMed = camis.CamisModel.fromRoots(isoRoots,isoRoots,isoRoots,isoRoots, gradient_threshold)
rIsoAsc = camis.CamisModel.fromRoots(CaRoots,CaRoots,CaRoots,CaRoots, gradient_threshold)

#AniCoLUT = camis.computeAniCoLUT(CdRoots, CaRoots, Cl1Roots, Cl2Roots, gradient_threshold)

"04 - Representation of CAMIS together with raw data"
camis.showCAMIS(CdRoots, CaRoots, Cl1Roots, Cl2Roots,beta,gradient,cost)
r1.showCAMIS()
rIsoMed.showCAMIS()
rIsoAsc.showCAMIS()

"05 - Saving CAMIS"
if fittingMode == 'DEM':
    camis.saveCamis('cuadriga_camis_dem.csv',CdRoots,CaRoots,Cl1Roots,Cl2Roots,\
                r1.anicoLUT)
if fittingMode == 'IMU':
    camis.saveCamis('cuadriga_camis_imu.csv',CdRoots,CaRoots,Cl1Roots,Cl2Roots,\
                r1.anicoLUT)
    camis.saveCamis('cuadriga_camis_imu_iso_med.csv',isoRoots,isoRoots,isoRoots,isoRoots,\
                rIsoMed.anicoLUT)
    camis.saveCamis('cuadriga_camis_imu_iso_asc.csv',CaRoots,CaRoots,CaRoots,CaRoots,\
                rIsoAsc.anicoLUT)




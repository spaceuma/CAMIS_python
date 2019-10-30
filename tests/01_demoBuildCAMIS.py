# -*- coding: utf-8 -*-
"""
DEMO script: building CAMIS based on Cuadriga experimental data
@author: J.Ricardo Sanchez Ibanez (ricardosan@uma.es), Github: @Ryk-San

This script is a reference demo to understand how CAMIS is built using data
obtained by a mobile robot.

"""

import numpy as np

import lib.cuadriga_reader as cr
import lib.camislib as camis
import lib.cost_mapping as costmapping
import matplotlib.pyplot as plt
import matplotlib.cm as cm


#hiRes_elevationMap = np.loadtxt(open("terrainData/UMARescueArea_10cmDEM.csv",\
#                                     "rb"), delimiter=" ", skiprows=0)
#hiRes_posX = np.loadtxt(open("terrainData/UMARescueArea_10cmPosX.csv",\
#                                     "rb"), delimiter=" ", skiprows=0)
#hiRes_posY = np.loadtxt(open("terrainData/UMARescueArea_10cmPosY.csv",\
#                                     "rb"), delimiter=" ", skiprows=0)
#
#hiRes = hiRes_posX[0,1] - hiRes_posX[0,0]
#env = costmapping.PDEM(hiRes_elevationMap, hiRes, .5, (hiRes_posX[0,0],hiRes_posY[0,0]))
#env.smoothMap(1.0)




#fig, axes = plt.subplots(constrained_layout=True)
#pathX, pathY, costPath, curvature, dHeading, segmentPath, dT = cr.getCuadrigaPath('cuadrigaData/201811/slopeTest09.txt')
#
#colors = cm.rainbow(np.linspace(0,1,len(dHeading)))
#axes.scatter(pathX, pathY, c = costPath, cmap = cm.rainbow)
##
#fig, axes = plt.subplots(constrained_layout=True)
#axes.plot(gradient)
###

fig, axes = plt.subplots(constrained_layout=True)
env.showMap('elevation',fig,axes)
#colors = cm.rainbow(np.linspace(0,1,len(dHeading)))
#axes.scatter(pathX-hiRes_posX[0,0], pathY-hiRes_posY[0,0], c = costPath, cmap = cm.rainbow)



"01 - The CSV files collected during preliminar experiments are loaded"
csvFiles = [
#            'cuadrigaData/20190531/JornadasRescate01.txt',
#            'cuadrigaData/20190531/JornadasRescate02.txt',
#            'cuadrigaData/20190531/JornadasRescate03.txt',
# #           'cuadrigaData/201811/slopeTest04.txt',
#            'cuadrigaData/201811/slopeTest05.txt',
#            'cuadrigaData/201811/slopeTest09.txt',
#            'cuadrigaData/201811/slopeTest10.txt',
            'cuadrigaData/20190624/2019_06_24_11_17_55.txt',
            'cuadrigaData/20190624/2019_06_24_11_22_21.txt',
            'cuadrigaData/20190624/2019_06_24_11_32_17.txt',
            'cuadrigaData/20190624/2019_06_24_11_36_06.txt',
            'cuadrigaData/20190624/2019_06_24_11_45_59.txt',
            'cuadrigaData/20190624/2019_06_24_11_49_29.txt',
            'cuadrigaData/20190624/2019_06_24_11_52_47.txt',
            'cuadrigaData/20190624/2019_06_24_11_57_35.txt'
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
#    G, Bx, By = camis.rpy2ab(Roll,Pitch,Yaw)
#    betaX = betaX + Bx
#    betaY = betaY + By
#    B = np.arctan2(By,Bx)
#    beta = beta + B.tolist()
#    currentCost = np.divide(Current,Speed)
#    cost = cost + currentCost.tolist()
#    gradient = gradient + G
    
    G, B = env.getBeta(posX - hiRes_posX[0,0], posY - hiRes_posY[0,0], heading)
    betaX = np.cos(B)
    betaY = np.sin(B)
    B = np.arctan2(betaY,betaX)
    beta = beta + B.tolist()
    currentCost = np.divide(Current,Speed)
    cost = cost + currentCost.tolist()
    G = np.abs(G)*180/np.pi
    gradient = gradient + G.tolist()
    
    
    

#beta.append(np.pi)
#cost.append(100)
#gradient.append(30)
#
#beta.append(0)
#cost.append(1000)
#gradient.append(30)
#
#beta.append(np.pi/2)
#cost.append(100)
#gradient.append(30)
#
#beta.append(-np.pi)
#cost.append(100)
#gradient.append(30)

gradient_threshold = 30.0

"03 - The functions Cd, Ca, Cl1 and Cl2 are obtained through a fitting process"
beta = [0 if np.isnan(x) else x for x in beta]
CdRoots, CaRoots, Cl1Roots, Cl2Roots = \
camis.computeDirCosts(gradient,beta,cost)

AniCoLUT = camis.computeAniCoLUT(CdRoots, CaRoots, Cl1Roots, Cl2Roots, 30.0)

"04 - Representation of CAMIS together with raw data"
camis.showCAMIS(CdRoots, CaRoots, Cl1Roots, Cl2Roots,beta,gradient,cost)

"05 - Saving CAMIS"
camis.saveCamis('cuadriga_camis.csv',CdRoots,CaRoots,Cl1Roots,Cl2Roots,\
                AniCoLUT)





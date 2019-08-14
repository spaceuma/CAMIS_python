# -*- coding: utf-8 -*-
"""
DEMO script: building CAMIS based on Cuadriga experimental data
@author: J.Ricardo Sanchez Ibanez (ricardosan@uma.es), Github: @Ryk-San

This script is a reference demo to understand how CAMIS is built using data
obtained by a mobile robot.

"""

import numpy as np

import lib.cuadriga_reader as cr
import lib.camis as camis


"01 - The CSV files collected during preliminar experiments are loaded"
csvFiles = [
            'cuadrigaData/20190531/JornadasRescate01.txt',
            'cuadrigaData/20190531/JornadasRescate02.txt',
            'cuadrigaData/20190531/JornadasRescate03.txt',
#            'cuadrigaData/201811/slopeTest04.txt',
            'cuadrigaData/201811/slopeTest05.txt',
            'cuadrigaData/201811/slopeTest09.txt',
            'cuadrigaData/201811/slopeTest10.txt',
            'cuadrigaData/20190624/2019_06_24_11_17_55.txt',
            'cuadrigaData/20190624/2019_06_24_11_22_21.txt',
            'cuadrigaData/20190624/2019_06_24_11_32_17.txt',
            'cuadrigaData/20190624/2019_06_24_11_36_06.txt',
            'cuadrigaData/20190624/2019_06_24_11_45_59.txt',
            'cuadrigaData/20190624/2019_06_24_11_49_29.txt',
            'cuadrigaData/20190624/2019_06_24_11_52_47.txt',
            'cuadrigaData/20190624/2019_06_24_11_57_35.txt'
            ]


"02 - CAMIS input data is extracted"
betaX = []
betaY = []
beta = []
cost = []
gradient = []
for file in csvFiles:
    " Cost C according to robot orientation Roll-Pitch-Yaw "
    Roll, Pitch, Yaw, Current, Speed = cr.readCuadrigaData(file)
    " Conversion from Roll-Pitch-Yaw to Gradient-Beta_angle "
    G, Bx, By = camis.rpy2ab(Roll,Pitch,Yaw)
    betaX = betaX + Bx
    betaY = betaY + By
    B = np.arctan2(By,Bx)
    beta = beta + B.tolist()
    currentCost = np.divide(Current,Speed)
    cost = cost + currentCost.tolist()
    gradient = gradient + G
    

beta.append(np.pi)
cost.append(100)
gradient.append(30)

beta.append(0)
cost.append(100)
gradient.append(30)

beta.append(np.pi/2)
cost.append(100)
gradient.append(30)

beta.append(-np.pi)
cost.append(100)
gradient.append(30)

"03 - The functions Cd, Ca, Cl1 and Cl2 are obtained through a fitting process"
beta = [0 if np.isnan(x) else x for x in beta]
CdRoots, CaRoots, Cl1Roots, Cl2Roots, AniCoLUT = \
camis.computeDirCosts(gradient,beta,cost)

"04 - Representation of CAMIS together with raw data"
camis.showCAMIS(CdRoots, CaRoots, Cl1Roots, Cl2Roots,beta,gradient,cost)

"05 - Saving CAMIS"
camis.saveCamis('cuadriga_camis.csv',CdRoots,CaRoots,Cl1Roots,Cl2Roots,\
                AniCoLUT)





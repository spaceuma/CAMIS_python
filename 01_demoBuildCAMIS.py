# -*- coding: utf-8 -*-
"""
DEMO script: building CAMIS based on Cuadriga experimental data
@author: J.Ricardo Sanchez Ibanez (ricardosan@uma.es), Github: @Ryk-San

This script is a reference demo to understand how CAMIS is built using data
obtained by a mobile robot.

"""

import numpy as np
from scipy.optimize import curve_fit

import lib.cuadriga_reader as cr
import lib.camis as camis
import csv


"01 - The CSV files collected during preliminar experiments are loaded"
csvFiles = ['cuadrigaData/20190531/JornadasRescate01.txt',
            'cuadrigaData/20190531/JornadasRescate02.txt',
            'cuadrigaData/20190531/JornadasRescate03.txt',
            'cuadrigaData/201811/slopeTest01.txt',
            'cuadrigaData/201811/slopeTest02.txt',
            'cuadrigaData/201811/slopeTest04.txt',
            'cuadrigaData/201811/slopeTest05.txt',
            'cuadrigaData/201811/slopeTest06.txt',
            'cuadrigaData/201811/slopeTest07.txt',
            'cuadrigaData/201811/slopeTest08.txt',
            'cuadrigaData/201811/slopeTest09.txt',
            'cuadrigaData/201811/slopeTest10.txt'
            ]


"02 - CAMIS input data is extracted"
betaX = []
betaY = []
beta = []
cost = []
gradient = []
for file in csvFiles:
    " Cost C according to robot orientation Roll-Pitch-Yaw "
    Roll, Pitch, Yaw, C = cr.readCuadrigaData(file)
    " Conversion from Roll-Pitch-Yaw to Gradient-Beta_angle "
    G, Bx, By = camis.rpy2ab(Roll,Pitch,Yaw)
    betaX = betaX + Bx
    betaY = betaY + By
    B = np.arctan2(By,Bx)
    beta = beta + B.tolist()
    cost = cost + C
    gradient = gradient + G
    
    
"03 - The functions Cd, Ca, Cl1 and Cl2 are obtained through a fitting process"
beta = [0 if np.isnan(x) else x for x in beta]
popt, pcov = curve_fit(camis.fittingCAMIS, (gradient,beta), cost)

"04 - Representation of CAMIS together with raw data"
camis.showCAMIS(popt,beta,gradient, cost)

"05 - Saving CAMIS"
with open('cuadriga_camis.csv', mode='w') as file:
    camis_writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    camis_writer.writerow(popt)





# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 15:23:37 2019

@author: rsanchez
"""

import numpy as np
import matplotlib.pyplot as plt
from context import camis
import copy

hiRes_elevationMap = np.loadtxt(\
                        open("data/terrainData/UMARescueArea_10cmDEM.csv",\
                             "rb"), delimiter=" ", skiprows=0)
hiRes_posX = np.loadtxt(open("data/terrainData/UMARescueArea_10cmPosX.csv",\
                                     "rb"), delimiter=" ", skiprows=0)
hiRes_posY = np.loadtxt(open("data/terrainData/UMARescueArea_10cmPosY.csv",\
                                     "rb"), delimiter=" ", skiprows=0)
hiRes = hiRes_posX[0,1] - hiRes_posX[0,0]


# =============================================================================
## CREATION OF CUADRIGA CAMIS
# =============================================================================
#exec(open("test_fittingCAMIS.py").read()) #Use twice, using DEM and IMU
demModel = camis.CamisModel.fromFile('cuadriga_camis_imu_iso_asc.csv')
imuModel = camis.CamisModel.fromFile('cuadriga_camis_imu.csv') 

envGo = camis.AnisotropicMap(hiRes_elevationMap[0:600,800:1200], hiRes, 0.5, (hiRes_posX[0,0],hiRes_posY[0,0]))
posA = np.asarray([30,50])
posB = np.asarray([15,5])

#envGo = camis.AnisotropicMap(hiRes_elevationMap[600:1100,0:600], hiRes, 0.5, (hiRes_posX[0,0],hiRes_posY[0,0]))
#posA = np.asarray([10,40])
#posB = np.asarray([50,45])

#envGo = camis.AnisotropicMap(hiRes_elevationMap[450:1250,700:1100], hiRes, 0.5, (hiRes_posX[0,0],hiRes_posY[0,0]))
#posA = np.asarray([5,50])
#posB = np.asarray([30,40])

envGo.smoothMap(2.0)
envGo_dem = copy.deepcopy(envGo)
envGo_imu = copy.deepcopy(envGo)

fig, axes = plt.subplots(constrained_layout=True)
envGo.showMap('elevation',fig,axes)


envGo_dem.computeVecCostMap(demModel)
envGo_imu.computeVecCostMap(imuModel)
envBack_dem = copy.deepcopy(envGo_dem)
envBack_imu = copy.deepcopy(envGo_imu)

envGo_dem.executePlanning(posB,posA)
envBack_dem.executePlanning(posA,posB)
envGo_imu.executePlanning(posB,posA)
envBack_imu.executePlanning(posA,posB)

fig, axes = plt.subplots(constrained_layout=True)
envGo.showMap('slope-deg',fig,axes)
envGo_dem.showPath(fig,axes,'r')
envBack_dem.showPath(fig,axes,'b')
envGo_imu.showPath(fig,axes,'m')
envBack_imu.showPath(fig,axes,'g')
axes.legend(('Isotropic-Go', 'Isotropic-Back', 'Anisotropic-Go', 'Anisotropic-Back'))

fig, axes = plt.subplots(constrained_layout=True)
envGo_dem.showPathData('cost',fig,axes,'r')
envBack_dem.showPathData('cost',fig,axes,'b')
envGo_imu.showPathData('cost',fig,axes,'m')
envBack_imu.showPathData('cost',fig,axes,'g')
axes.legend(('Isotropic-Go', 'Isotropic-Back', 'Anisotropic-Go', 'Anisotropic-Back'))

fig, axes = plt.subplots(constrained_layout=True)
envGo_dem.showPathData('total-cost',fig,axes,'r')
envGo_dem.showPathData('total-cost-estimated',fig,axes,'r')
envBack_dem.showPathData('total-cost',fig,axes,'b')
envBack_dem.showPathData('total-cost-estimated',fig,axes,'b')
envGo_imu.showPathData('total-cost',fig,axes,'m')
envGo_imu.showPathData('total-cost-estimated',fig,axes,'m')
envBack_imu.showPathData('total-cost',fig,axes,'g')
envBack_imu.showPathData('total-cost-estimated',fig,axes,'g')
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
demModel = camis.CamisModel.fromFile('cuadriga_camis_dem_iso_asc.csv')
imuModel = camis.CamisModel.fromFile('cuadriga_camis_dem.csv') 

#envGo = camis.AnisotropicMap(hiRes_elevationMap[0:550,800:1200], hiRes, 0.5, (hiRes_posX[0,0],hiRes_posY[0,0]))
#posA = np.asarray([30,50])
#posB = np.asarray([15,5])

envGo = camis.AnisotropicMap(hiRes_elevationMap[600:1100,0:600], hiRes, 0.5, (hiRes_posX[0,0],hiRes_posY[0,0]))
posA = np.asarray([10,40])
posB = np.asarray([50,45])
#
#envGo = camis.AnisotropicMap(hiRes_elevationMap[450:1250,700:1100], hiRes, 0.5, (hiRes_posX[0,0],hiRes_posY[0,0]))
#posA = np.asarray([5,50])
#posB = np.asarray([30,40])

envGo.smoothMap(2.0)
envGo.processSlope(10.0)

envGo_dem = copy.deepcopy(envGo)
envGo_imu = copy.deepcopy(envGo)

fig, axes = plt.subplots(constrained_layout=True)
envGo.showMap('elevation',fig,axes)


envGo_dem.computeVecCostMap(demModel)
envGo_imu.computeVecCostMap(imuModel)
envBack_dem = copy.deepcopy(envGo_dem)
envGo_bidem = copy.deepcopy(envGo_dem)
envBack_bidem = copy.deepcopy(envGo_dem)
envBack_imu = copy.deepcopy(envGo_imu)
envGo_biimu = copy.deepcopy(envGo_imu)
envBack_biimu = copy.deepcopy(envGo_imu)

envGo_bidem.executeBiPlanning(posB,posA)
envBack_bidem.executeBiPlanning(posA,posB)
envGo_biimu.executeBiPlanning(posB,posA)
envBack_biimu.executeBiPlanning(posA,posB)

envGo_dem.executePlanning(posB,posA)
envBack_dem.executePlanning(posA,posB)
envGo_imu.executePlanning(posB,posA)
envBack_imu.executePlanning(posA,posB)

fig, axes = plt.subplots(constrained_layout=True)
envGo.showMap('old-elevation',fig,axes)
envGo_bidem.showPath(fig,axes,'r','solid')
envGo_dem.showPath(fig,axes,'r','dotted')
envBack_bidem.showPath(fig,axes,'b','solid')
envBack_dem.showPath(fig,axes,'b','dotted')
envGo_biimu.showPath(fig,axes,'m','solid')
envGo_imu.showPath(fig,axes,'m','dotted')
envBack_biimu.showPath(fig,axes,'g','solid')
envBack_imu.showPath(fig,axes,'g','dotted')
axes.legend(('Isotropic-Go', 'Isotropic-Back', 'Anisotropic-Go', 'Anisotropic-Back'))

fig, axes = plt.subplots(constrained_layout=True)
envGo_dem.showPathData('cost',fig,axes,'r')
envBack_dem.showPathData('cost',fig,axes,'b')
envGo_imu.showPathData('cost',fig,axes,'m')
envBack_imu.showPathData('cost',fig,axes,'g')
axes.legend(('Isotropic-Go', 'Isotropic-Back', 'Anisotropic-Go', 'Anisotropic-Back'))

fig, axes = plt.subplots(constrained_layout=True)
envGo_bidem.showPathData('total-cost-estimated',fig,axes,'r')
envGo_dem.showPathData('total-cost-estimated',fig,axes,'r')
envBack_bidem.showPathData('total-cost-estimated',fig,axes,'b')
envBack_dem.showPathData('total-cost-estimated',fig,axes,'b')
envGo_biimu.showPathData('total-cost-estimated',fig,axes,'m')
envGo_imu.showPathData('total-cost-estimated',fig,axes,'m')
envBack_biimu.showPathData('total-cost-estimated',fig,axes,'g')
envBack_imu.showPathData('total-cost-estimated',fig,axes,'g')

#envGo.show3dDEM()
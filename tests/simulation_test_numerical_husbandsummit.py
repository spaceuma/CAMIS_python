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
# hiRes_elevationMap = np.loadtxt(\
#                         open("data/terrainData/husbandSummit/husbandsummit_elevationMap.txt",\
#                         "rb"), delimiter=" ", skiprows=0)
# hiRes_elevationMap = hiRes_elevationMap[0:64,:] 
# hiRes_elevationMap = np.loadtxt(\
#                           open("data/terrainData/victoriaCrater/victoriacrater_elevationMap.txt",\
#                           "rb"), delimiter=" ", skiprows=0)
# hiRes_elevationMap = hiRes_elevationMap[170:250,0:80]

   
hiRes_elevationMap = np.loadtxt(\
                          open("data/terrainData/crater/crater_elevationMap.txt",\
                          "rb"), delimiter=" ", skiprows=0)
# hiRes_elevationMap = hiRes_elevationMap[0:32,:] 

hiRes_elevationMap = hiRes_elevationMap - hiRes_elevationMap.min()
hiRes_elevationMap = hiRes_elevationMap/2.0



hiRes = 1.0
offset = [0,0]



hexRes = 0.5

occupancy_radius = 0.5
tracking_error = 0.5
env = camis.AnisotropicMap(hiRes_elevationMap, hiRes, hexRes,\
                               offset, occupancy_radius, tracking_error)

posA = np.asarray([10,10])
posB = np.asarray([55,50])
# posB = np.asarray([58,40])

plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'

fig, ax = plt.subplots(constrained_layout=True)
cc = ax.contourf(hiRes_elevationMap, 
                     50, cmap = cm.gist_earth, extend='both')
ax.contour(hiRes_elevationMap, 
               10, colors = 'k', alpha=.3)
    

env.showHexElevationMap(16)

env.showSqSlopeGradientMap(16)

def getMapLists(camisInput):
    iso_model = copy.deepcopy(camisInput)
    iso_model.setAsIsotropic()
    env_aniso = [None] * 2
    env_iso = [None] * 2
    env.computeVecCostMap(camisInput)
    enviso = copy.deepcopy(env)
    enviso.computeVecCostMap(iso_model)
    for i in range(2):
        env_aniso[i] = copy.deepcopy(env)
        env_iso[i] = copy.deepcopy(enviso)
    return env_aniso, env_iso

def computeAllPlannings(anisoMapList):
    anisoMapList[0].executeHexBiPlanning(posB,posA, nbUpdate = True, 
                                                  anisoSearch = 'single',
                                                  dirPolicy = 'bidir')
    anisoMapList[1].executeHexBiPlanning(posA,posB, nbUpdate = True, 
                                                  anisoSearch = 'single',
                                                  dirPolicy = 'bidir')



with open("data/thesis/cuadriga_aniso_01.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_01 = camis.CamisDrivingModel(cuadriga_data)
aniso_01.showCAMIS(20)
env_CUAD01_scene01, env_isoCUAD01_scene01 = getMapLists(aniso_01)
computeAllPlannings(env_CUAD01_scene01)
computeAllPlannings(env_isoCUAD01_scene01)

with open("data/thesis/cuadriga_aniso_02.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_02 = camis.CamisDrivingModel(cuadriga_data)
aniso_02.showCAMIS(20)
env_CUAD02_scene01, env_isoCUAD02_scene01 = getMapLists(aniso_02)
computeAllPlannings(env_CUAD02_scene01)
computeAllPlannings(env_isoCUAD02_scene01)

with open("data/thesis/cuadriga_aniso_03.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_03 = camis.CamisDrivingModel(cuadriga_data)
aniso_03.showCAMIS(20)
env_CUAD03_scene01, env_isoCUAD03_scene01 = getMapLists(aniso_03)
computeAllPlannings(env_CUAD03_scene01)
computeAllPlannings(env_isoCUAD03_scene01)

with open("data/thesis/cuadriga_aniso_04.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_04 = camis.CamisDrivingModel(cuadriga_data)
aniso_04.showCAMIS(20)
env_CUAD04_scene01, env_isoCUAD04_scene01 = getMapLists(aniso_04)
computeAllPlannings(env_CUAD04_scene01)
computeAllPlannings(env_isoCUAD04_scene01)

with open("data/thesis/cuadriga_aniso_05.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_05 = camis.CamisDrivingModel(cuadriga_data)
# aniso_05.showCAMIS(25)
env_CUAD05_scene01, env_isoCUAD05_scene01 = getMapLists(aniso_05)
computeAllPlannings(env_CUAD05_scene01)
computeAllPlannings(env_isoCUAD05_scene01)

with open("data/thesis/cuadriga_aniso_06.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_06 = camis.CamisDrivingModel(cuadriga_data)
aniso_06.showCAMIS(20)
env_CUAD06_scene01, env_isoCUAD06_scene01 = getMapLists(aniso_06)
computeAllPlannings(env_CUAD06_scene01)
computeAllPlannings(env_isoCUAD06_scene01)

with open("data/thesis/cuadriga_aniso_07.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_07 = camis.CamisDrivingModel(cuadriga_data)
aniso_07.showCAMIS(20)
env_CUAD07_scene01, env_isoCUAD07_scene01 = getMapLists(aniso_07)
computeAllPlannings(env_CUAD07_scene01)
computeAllPlannings(env_isoCUAD07_scene01)

with open("data/thesis/cuadriga_aniso_08.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_08 = camis.CamisDrivingModel(cuadriga_data)
# aniso_08.showCAMIS(25)
env_CUAD08_scene01, env_isoCUAD08_scene01 = getMapLists(aniso_08)
computeAllPlannings(env_CUAD08_scene01)
computeAllPlannings(env_isoCUAD08_scene01)

with open("data/thesis/cuadriga_aniso_09.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_09 = camis.CamisDrivingModel(cuadriga_data)
# aniso_09.showCAMIS(25)
env_CUAD09_scene01, env_isoCUAD09_scene01 = getMapLists(aniso_09)
computeAllPlannings(env_CUAD09_scene01)
computeAllPlannings(env_isoCUAD09_scene01)


with open("data/thesis/cuadriga_aniso_10.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_10 = camis.CamisDrivingModel(cuadriga_data)
# aniso_10.showCAMIS(25)
env_CUAD10_scene01, env_isoCUAD10_scene01 = getMapLists(aniso_10)
computeAllPlannings(env_CUAD10_scene01)
computeAllPlannings(env_isoCUAD10_scene01)

with open("data/thesis/cuadriga_aniso_11.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_11 = camis.CamisDrivingModel(cuadriga_data)
# aniso_11.showCAMIS(25)
env_CUAD11_scene01, env_isoCUAD11_scene01 = getMapLists(aniso_11)
computeAllPlannings(env_CUAD11_scene01)
computeAllPlannings(env_isoCUAD11_scene01)

with open("data/thesis/cuadriga_aniso_12.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_12 = camis.CamisDrivingModel(cuadriga_data)
# aniso_12.showCAMIS(25)
env_CUAD12_scene01, env_isoCUAD12_scene01 = getMapLists(aniso_12)
computeAllPlannings(env_CUAD12_scene01)
computeAllPlannings(env_isoCUAD12_scene01)

# env_CUAD01_scene01[0].executeHexBiPlanning(posB, posA, nbUpdate = True, 
#                                            anisoSearch = 'single',
#                                            dirPolicy = 'bidir')
# env_CUAD01_scene01[0].showHexBiTmaps()

# env_CUAD01_scene01[1].executeHexBiPlanning(posA, posB, nbUpdate = True, 
#                                            anisoSearch = 'single',
#                                            dirPolicy = 'bidir')
# env_CUAD01_scene01[1].showHexBiTmaps()


# env_isoCUAD01_scene01[0].executeHexBiPlanning(posB, posA, nbUpdate = True, 
#                                              anisoSearch = 'single',
#                                              dirPolicy = 'bidir')
# env_isoCUAD01_scene01[0].showHexBiTmaps()

# env_isoCUAD01_scene01[1].executeHexBiPlanning(posA, posB, nbUpdate = True, 
#                                              anisoSearch = 'single',
#                                              dirPolicy = 'bidir')
# env_isoCUAD01_scene01[1].showHexBiTmaps()


##### 3D PLOT

from matplotlib.colors import LightSource

XX = np.zeros_like(hiRes_elevationMap)
YY = np.zeros_like(hiRes_elevationMap)
for i in range(hiRes_elevationMap.shape[0]):
    for j in range(hiRes_elevationMap.shape[1]):
        XX[j,i] = i
        YY[j,i] = j




plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
# fig = plt.figure(figsize=(14.0, 7.0))
ax1 = plt.figure(figsize=(7.5, 6.5)).add_subplot(projection='3d')
ax1.view_init(elev=80.0, azim= 125.0)


ls = LightSource(270, 5)
rgb = ls.shade(hiRes_elevationMap, cmap=cm.gist_earth, vert_exag=1.0, blend_mode='soft')
surf1 = ax1.plot_surface(XX,YY,hiRes_elevationMap, rstride=1, cstride=1, facecolors=rgb,\
                        linewidth=0, antialiased=True, shade=False,\
                        vmin = 0.0, vmax = 1.0, rasterized=True)

anisoresults = [env_CUAD07_scene01[0], env_CUAD01_scene01[0],
                env_CUAD08_scene01[0], env_CUAD02_scene01[0],
                env_CUAD09_scene01[0], env_CUAD03_scene01[0],
                env_CUAD10_scene01[0], env_CUAD04_scene01[0],
                env_CUAD11_scene01[0], env_CUAD05_scene01[0],
                env_CUAD12_scene01[0], env_CUAD06_scene01[0]]

isoresults = [env_isoCUAD07_scene01[0], env_isoCUAD01_scene01[0],
              env_isoCUAD08_scene01[0], env_isoCUAD02_scene01[0],
              env_isoCUAD09_scene01[0], env_isoCUAD03_scene01[0],
              env_isoCUAD10_scene01[0], env_isoCUAD04_scene01[0],
              env_isoCUAD11_scene01[0], env_isoCUAD05_scene01[0],
              env_isoCUAD12_scene01[0], env_isoCUAD06_scene01[0]]

for i,scene in enumerate(anisoresults):
    ax1.plot(scene.path[:,0], scene.path[:,1], scene.pathElevation,
             linestyle='solid',color = plt.cm.jet(float(i)/12.0))

for i,scene in enumerate(isoresults):
    ax1.plot(scene.path[:,0], scene.path[:,1], scene.pathElevation,
             linestyle='dashed',color = plt.cm.jet(float(i)/12.0))   

ax1.plot(env_isoCUAD07_scene01[0].path[:,0],
         env_isoCUAD07_scene01[0].path[:,1], 
         env_isoCUAD07_scene01[0].pathElevation, linestyle='dotted',
         color = plt.cm.jet(1.0/12.0))  
ax1.plot(env_isoCUAD01_scene01[0].path[:,0], 
         env_isoCUAD01_scene01[0].path[:,1],
         env_isoCUAD01_scene01[0].pathElevation, linestyle='dotted',
         color = plt.cm.jet(2.0/12.0))    
ax1.plot(env_isoCUAD08_scene01[0].path[:,0],
         env_isoCUAD08_scene01[0].path[:,1], 
         env_isoCUAD08_scene01[0].pathElevation, linestyle='dotted',
         color = plt.cm.jet(3.0/12.0))
ax1.plot(env_isoCUAD02_scene01[0].path[:,0], 
         env_isoCUAD02_scene01[0].path[:,1], \
         env_isoCUAD02_scene01[0].pathElevation,linestyle='dotted',
         color = plt.cm.jet(4.0/12.0))
ax1.plot(env_isoCUAD09_scene01[0].path[:,0], 
         env_isoCUAD09_scene01[0].path[:,1], \
         env_isoCUAD09_scene01[0].pathElevation,linestyle='dotted',
         color = plt.cm.jet(5.0/12.0))
ax1.plot(env_isoCUAD03_scene01[0].path[:,0], env_isoCUAD03_scene01[0].path[:,1], \
        env_isoCUAD03_scene01[0].pathElevation,linestyle='dotted',color = 'orange')  
ax1.plot(env_isoCUAD04_scene01[0].path[:,0], env_isoCUAD04_scene01[0].path[:,1], \
        env_isoCUAD04_scene01[0].pathElevation,linestyle='dotted',color = 'b')
ax1.plot(env_isoCUAD05_scene01[0].path[:,0], env_isoCUAD05_scene01[0].path[:,1], \
        env_isoCUAD05_scene01[0].pathElevation,linestyle='dotted',color = 'c')
ax1.plot(env_isoCUAD06_scene01[0].path[:,0], env_isoCUAD06_scene01[0].path[:,1], \
        env_isoCUAD06_scene01[0].pathElevation,linestyle='dotted',color = 'lime')

    
ax1.set_xlabel('X-axis [m]',labelpad=10,fontsize='x-large')
ax1.set_ylabel('Y-axis [m]',labelpad=10,fontsize='x-large')
ax1.set_zlabel('Z-axis [m]',labelpad=0,fontsize='x-large')
ax1.tick_params(axis="x",direction="in", pad=0, labelsize = 'x-large')
ax1.tick_params(axis="y",direction="in", pad=0, labelsize = 'x-large')
ax1.tick_params(axis="z",direction="in", pad=0, labelsize = 'x-large')
ax1.set_zticks([0.0, 1.0, 2.0])
# ax1.view_init(78.0, -150.0)
# ax1.set_xlim([0,50])
# ax1.set_ylim([0,50])
ax1.set_zlim([0,2.0])
ax1.set_facecolor('w')
plt.subplots_adjust(left = 0.01, right = 0.9, bottom = 0.0, top = 1.0, \
                        wspace = 0.0, hspace = 0.0)
fig.tight_layout()


ax2 = plt.figure(figsize=(7.0, 7.0)).add_subplot(projection='3d')
ax2.view_init(elev=80.0, azim= 125.0)
surf2 = ax2.plot_surface(XX,YY,hiRes_elevationMap, rstride=1, cstride=1, facecolors=rgb,\
                        linewidth=0, antialiased=True, shade=False,\
                        vmin = 0.0, vmax = 1.0, rasterized=True)
    
ax2.plot(env_CUAD01_scene01[1].path[:,0], env_CUAD01_scene01[1].path[:,1], \
        env_CUAD01_scene01[1].pathElevation,linestyle='solid',color = 'r')
ax2.plot(env_CUAD02_scene01[1].path[:,0], env_CUAD02_scene01[1].path[:,1], \
        env_CUAD02_scene01[1].pathElevation,linestyle='solid',color = 'm')
ax2.plot(env_CUAD03_scene01[1].path[:,0], env_CUAD03_scene01[1].path[:,1], \
        env_CUAD03_scene01[1].pathElevation,linestyle='solid',color = 'orange')  
ax2.plot(env_CUAD04_scene01[1].path[:,0], env_CUAD04_scene01[1].path[:,1], \
        env_CUAD04_scene01[1].pathElevation,linestyle='solid',color = 'b')
ax2.plot(env_CUAD05_scene01[1].path[:,0], env_CUAD05_scene01[1].path[:,1], \
        env_CUAD05_scene01[1].pathElevation,linestyle='solid',color = 'c')
ax2.plot(env_CUAD06_scene01[1].path[:,0], env_CUAD06_scene01[1].path[:,1], \
        env_CUAD06_scene01[1].pathElevation,linestyle='solid',color = 'lime')    
ax2.plot(env_isoCUAD01_scene01[0].path[:,0], env_isoCUAD01_scene01[0].path[:,1], \
        env_isoCUAD01_scene01[0].pathElevation,linestyle='dotted',color = 'r')    
ax2.plot(env_isoCUAD02_scene01[0].path[:,0], env_isoCUAD02_scene01[0].path[:,1], \
        env_isoCUAD02_scene01[0].pathElevation,linestyle='dotted',color = 'm')
ax2.plot(env_isoCUAD03_scene01[0].path[:,0], env_isoCUAD03_scene01[0].path[:,1], \
        env_isoCUAD03_scene01[0].pathElevation,linestyle='dotted',color = 'orange')  
ax2.plot(env_isoCUAD04_scene01[0].path[:,0], env_isoCUAD04_scene01[0].path[:,1], \
        env_isoCUAD04_scene01[0].pathElevation,linestyle='dotted',color = 'b')
ax2.plot(env_isoCUAD05_scene01[0].path[:,0], env_isoCUAD05_scene01[0].path[:,1], \
        env_isoCUAD05_scene01[0].pathElevation,linestyle='dotted',color = 'c')
ax2.plot(env_isoCUAD06_scene01[0].path[:,0], env_isoCUAD06_scene01[0].path[:,1], \
        env_isoCUAD06_scene01[0].pathElevation,linestyle='dotted',color = 'lime')
   
ax2.set_xlabel('X-axis [m]',labelpad=-8)
ax2.set_ylabel('Y-axis [m]',labelpad=-8)
ax2.set_zlabel('Z-axis [m]',labelpad=-6)
# ax2.tick_params(axis="x",direction="in", pad=-5)
# ax2.tick_params(axis="y",direction="in", pad=-5)
# ax2.tick_params(axis="z",direction="in", pad=-2)
# ax2.set_zticks([-0.6, 0.1, 0.8])
# ax2.view_init(78.0, -150.0)
# ax2.set_xlim([0,50])
# ax2.set_ylim([0,50])
ax2.set_zlim([0,2.0])
ax2.set_facecolor('w')


plt.subplots_adjust(left = 0.01, right = 1.0, bottom = 0.0, top = 1.0, \
                        wspace = 0.0, hspace = 0.0)

    
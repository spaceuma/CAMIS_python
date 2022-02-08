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
aniso_01.showCAMIS(20, iso = True)
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

fig = plt.figure(figsize=(8.5, 6.5))
ax1 = fig.add_subplot(projection='3d')
ax1.view_init(elev=80.0, azim= 135.0)

ls = LightSource(270, 5)
rgb = ls.shade(hiRes_elevationMap, cmap=cm.gist_earth, vert_exag=1.0, 
               blend_mode='soft')
surf1 = ax1.plot_surface(XX,YY,hiRes_elevationMap, rstride=1, cstride=1, 
                         facecolors=rgb, linewidth=0, antialiased=True, 
                         shade=False, vmin = 0.0, vmax = 1.0, rasterized=True)

anisowheelresults = [env_CUAD07_scene01[0], env_CUAD01_scene01[0],
                     env_CUAD08_scene01[0], env_CUAD02_scene01[0],
                     env_CUAD09_scene01[0], env_CUAD03_scene01[0]]

anisotrackresults = [env_CUAD10_scene01[0], env_CUAD04_scene01[0],
                     env_CUAD11_scene01[0], env_CUAD05_scene01[0],
                     env_CUAD12_scene01[0], env_CUAD06_scene01[0]]

isowheelresults = [env_isoCUAD07_scene01[0], env_isoCUAD01_scene01[0],
                   env_isoCUAD08_scene01[0], env_isoCUAD02_scene01[0],
                   env_isoCUAD09_scene01[0], env_isoCUAD03_scene01[0]]

isotrackresults = [env_isoCUAD10_scene01[0], env_isoCUAD04_scene01[0],
                   env_isoCUAD11_scene01[0], env_isoCUAD05_scene01[0],
                   env_isoCUAD12_scene01[0], env_isoCUAD06_scene01[0]]

for i,scene in enumerate(anisowheelresults):
    ax1.plot(scene.path[:,0], scene.path[:,1], scene.pathElevation,
             linestyle='solid',color = \
             plt.cm.jet(float(i)/float(len(anisowheelresults))))
        
for i,scene in enumerate(anisotrackresults):
    ax1.plot(scene.path[:,0], scene.path[:,1], scene.pathElevation,
             linestyle='dashdot',color = \
             plt.cm.jet(float(i)/float(len(anisotrackresults))))

for i,scene in enumerate(isowheelresults):
    ax1.plot(scene.path[:,0], scene.path[:,1], scene.pathElevation,
             linestyle='dotted',color = \
             plt.cm.jet(float(i)/float(len(isowheelresults))))    

for i,scene in enumerate(isotrackresults):
    ax1.plot(scene.path[:,0], scene.path[:,1], scene.pathElevation,
             linestyle='dashed',color = \
             plt.cm.jet(float(i)/float(len(isotrackresults)))) 
        
ax1.set_xlabel('X-axis [m]',labelpad=10,fontsize='x-large')
ax1.set_ylabel('Y-axis [m]',labelpad=10,fontsize='x-large')
ax1.set_zlabel('Z-axis [m]',labelpad=0,fontsize='x-large')
ax1.tick_params(axis="x",direction="in", pad=0, labelsize = 'x-large')
ax1.tick_params(axis="y",direction="in", pad=0, labelsize = 'x-large')
ax1.tick_params(axis="z",direction="in", pad=0, labelsize = 'x-large')
ax1.set_zticks([0.0, 1.0, 2.0])
ax1.set_zlim([0,2.0])
ax1.set_facecolor('w')
plt.subplots_adjust(left = 0.05, right = 0.9, bottom = 0.0, top = 1.0, \
                        wspace = 0.0, hspace = 0.0)

ax1.text(posA[0]-1,posA[1]-1,2.0,'$X_o$',color='k', size=15)    
ax1.text(posB[0]+4,posB[1]+4,2.0,'$X_g$',color='k', size=15)
  
# fig.tight_layout()

legend_ax_wheel_aniso = fig.add_axes([0.0, 0.8, 0.2, 0.2])
legend_ax_track_aniso = fig.add_axes([0.0, 0.0, 0.2, 0.45])
legend_ax_wheel_iso = fig.add_axes([0.8, 0.8, 0.2, 0.2])
legend_ax_track_iso = fig.add_axes([0.8, 0.0, 0.2, 0.45])

legend_ax_wheel_aniso.plot(0,0, alpha = 0, label = '$\mathbf{Anisotropic}$')
legend_ax_wheel_aniso.plot(0,0, alpha = 0, label = '$\mathbf{Wheel \ Model}$')
legend_ax_wheel_aniso.plot(0,0, alpha = 0, label = '$\mathbf{σ_{ij} = 0.07 e^{0.1 α_{ij}}}$')
legend_ax_track_aniso.plot(0,0, alpha = 0, label = '$\mathbf{Anisotropic}$')
legend_ax_track_aniso.plot(0,0, alpha = 0, label = '$\mathbf{Track \ Model}$')
legend_ax_track_aniso.plot(0,0, alpha = 0, label = '$\mathbf{σ_{ij} = 0.04 e^{0.07 α_{ij}}}$')
legend_ax_wheel_iso.plot(0,0, alpha = 0, label = '$\mathbf{Isotropic}$')
legend_ax_wheel_iso.plot(0,0, alpha = 0, label = '$\mathbf{Wheel \ Model}$')
legend_ax_wheel_iso.plot(0,0, alpha = 0, label = '$\mathbf{σ_{ij} = 0.07 e^{0.1 α_{ij}}}$')
legend_ax_track_iso.plot(0,0, alpha = 0, label = '$\mathbf{Isotropic}$')
legend_ax_track_iso.plot(0,0, alpha = 0, label = '$\mathbf{Track \ Model}$')
legend_ax_track_iso.plot(0,0, alpha = 0, label = '$\mathbf{σ_{ij} = 0.04 e^{0.07 α_{ij}}}$')


legends = [legend_ax_wheel_aniso, legend_ax_track_aniso, 
           legend_ax_wheel_iso, legend_ax_track_iso]

labels = ['$ρ_{ij} = 0.15$', '$ρ_{ij} = 0.3$', '$ρ_{ij} = 0.45$',
          '$ρ_{ij} = 0.6$', '$ρ_{ij} = 0.75$', '$ρ_{ij} = 0.9$']

styles = ['solid','dashdot','dotted','dashed']

for j,leg in enumerate(legends):
    for i,l in enumerate(labels):
        leg.plot(0,0,c=plt.cm.jet(float(i)/float(len(labels))), 
                 linestyle = styles[j], label=l)
    leg.grid(visible=False)
    leg.axis('off')
    leg.legend(fontsize=12,ncol=1)
    
    
    
### ENERGY 
    
plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
coeffsLabels = ['$ρ_{ij}= 0.15$','$ρ_{ij} = 0.3$',\
                '$ρ_{ij} = 0.45$','$ρ_{ij} = 0.6$',\
                '$ρ_{ij} = 0.75$','$ρ_{ij} = 0.9$']
anisoTotalCostWheel = [0,0,0,0,0,0]
anisoTotalCostEstimatedWheel = [0,0,0,0,0,0]
anisoNumUpdatesWheel = [0,0,0,0,0,0]

isoTotalCostWheel = [0,0,0,0,0,0]
isoTotalCostEstimatedWheel = [0,0,0,0,0,0]
isoNumUpdatesWheel = [0,0,0,0,0,0]


anisoTotalCostTrack = [0,0,0,0,0,0]
anisoTotalCostEstimatedTrack = [0,0,0,0,0,0]
isoTotalCostTrack = [0,0,0,0,0,0]
isoTotalCostEstimatedTrack = [0,0,0,0,0,0]
anisoNumUpdatesTrack = [0,0,0,0,0,0]
isoNumUpdatesTrack = [0,0,0,0,0,0]
anisoTR = [0,0]
isoTR = [0,0]

for i,res in enumerate(anisowheelresults):
    anisoTotalCostWheel[i] = res.pathComputedTotalCost[-1]
    anisoTotalCostEstimatedWheel[i] = res.pathEstimatedTotalCost[0]
    anisoNumUpdatesWheel[i] = res.numUpdates

for i,res in enumerate(isowheelresults):
    isoTotalCostWheel[i] =  res.pathComputedTotalCost[-1]
    isoTotalCostEstimatedWheel[i] = res.pathEstimatedTotalCost[0]
    isoNumUpdatesWheel[i] = res.numUpdates
    
for i,res in enumerate(anisotrackresults):
    anisoTotalCostTrack[i] = res.pathComputedTotalCost[-1]
    anisoTotalCostEstimatedTrack[i] = res.pathEstimatedTotalCost[0]
    anisoNumUpdatesTrack[i] = res.numUpdates
    
for i,res in enumerate(isotrackresults):
    isoTotalCostTrack[i] = res.pathComputedTotalCost[-1]
    isoTotalCostEstimatedTrack[i] = res.pathEstimatedTotalCost[0]
    isoNumUpdatesTrack[i] = res.numUpdates

x = np.arange(len(coeffsLabels))  # the label locations
x2 = np.arange(2)+1
width = 0.48  # the width of the bars

fig, ax = plt.subplots(figsize=(10,1.5), constrained_layout=True)
#rects3 = ax.bar(x2 - 0.45/2, anisoTR, 0.45, label='Isotropic (ρ = 0.8)', color='lime')
#rects4 = ax.bar(x2 + 0.45/2, isoTR, 0.45, label='Isotropic (ρ = 0.8)', color='g')
rects1 = ax.bar(x - width/2, anisoTotalCostTrack, width, label='Isotropic (ρ_{ij} = 0.8)', color='r')
rects2 = ax.bar(x + width/2, isoTotalCostTrack, width, label='Isotropic (ρ_{ij} = 0.8)', color = 'b')

gainWheel = []

for i,rect in enumerate(rects1):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T),
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, -2),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='top',color = 'w', fontsize = 12)
for i,rect in enumerate(rects2):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T),
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, -2),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='top',color = 'w', fontsize = 12)
for i,rect in enumerate(rects2):
        isoT = rect.get_height()
        anisoT = rects1[i].get_height()
        gain = (isoT - anisoT)/isoT * 100
        gainWheel.append(gain)
        ax.annotate('-{0:.2f}'.format(gain) + '%',
                    xy=(rect.get_x(), isoT),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize = 12) 
        
ax.grid(True, which='both')   
#autolabel(rects1)
#autolabel(rects2)
ax.set_ylabel('Total Cost [As]', fontsize = 12)
ax.set_ylim([0,1200])
ax.set_xlim([-0.75,5.75])
#ax.set_xlabel('CAMIS')
ax.set_xticks(x)
ax.set_xticklabels(coeffsLabels, fontsize = 12)
ax.legend(('Anisotropic','Isotropic'), fontsize = 12)
plt.minorticks_on()  
plt.show()

fig, ax = plt.subplots(figsize=(5,4), constrained_layout=True)
ax.plot(np.divide(anisoNumUpdatesWheel,isoNumUpdatesWheel))
ax.plot(np.divide(anisoNumUpdatesTrack,isoNumUpdatesTrack))
ax.set_ylabel('Number of total cost updates \n in OUM with respect to FMM', fontsize = 14)
ax.set_xlabel('Specific resistance $ρ_{ij}$', fontsize = 14)
ax.legend(('Wheel Model','Track Model'), fontsize = 14)

fig, ax = plt.subplots(figsize=(5,4), constrained_layout=True)
ax.plot(np.divide(np.subtract(anisoTotalCostWheel,isoTotalCostWheel),
                  anisoTotalCostWheel)*100)
ax.plot(np.divide(np.subtract(anisoTotalCostTrack,isoTotalCostTrack),
                  anisoTotalCostTrack)*100)
ax.set_ylabel('Reduction of total cost \n considering anisotropic cost \n [%]', fontsize = 14)
ax.set_xlabel('Specific resistance $ρ_{ij}$', fontsize = 14)
# ax.set_xticks(np.arange(len(coeffsLabels)))
# ax.set_xticklabels(coeffsLabels)
ax.legend(('Wheel Model','Track Model'), fontsize = 14)


### Total Cost Maps
env_CUAD07_scene01[0].showHexBiTmaps()
env_isoCUAD07_scene01[0].showHexBiTmaps()

env_CUAD01_scene01[0].showHexBiTmaps()
env_isoCUAD01_scene01[0].showHexBiTmaps()
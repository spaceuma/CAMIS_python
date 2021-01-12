# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:35:43 2020

@author: Richi
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from context import camis
import copy
import zipfile
import yaml
try:
    from scipy import signal
except:
    raise ImportError('ERROR: scipy module could not be imported')
import plotly.graph_objects as go

# =============================================================================
## LOADING DEM
# =============================================================================

#v = x.*exp(-x.^2-y.^2-z.^2);

xx = np.linspace(-4.5,4.5,50)
yy = np.linspace(-4.5,4.5,50)
xmesh, ymesh = np.meshgrid(xx,yy)

z = np.ones_like(xmesh)
for j,y in enumerate(yy):
    for i,x in enumerate(xx):
        z[j,i] =  3*(1-x)**2*np.exp(-(x**2) - (y+1)**2) - \
        10*(x/5 - x**3 - y**5)*np.exp(-x**2-y**2) - 1/3*np.exp(-(x+1)**2 - y**2) 

xmesh *= 10.0
ymesh *= 10.0
z = z/10.0

plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, ax1 = plt.subplots(figsize=(5, 4),nrows = 1, ncols = 1, constrained_layout=True)
cc = ax1.scatter(xmesh, ymesh, c = z, 
                 cmap=cm.gist_earth,s=16.0)
cbar = fig.colorbar(cc, ax=ax1,shrink=0.6)

hiRes = 1.0
hexRes = 0.5
offset = (0,0)
occupancy_radius = 0.5
tracking_error = 0.5
posA = np.asarray([20,45])
posB = np.asarray([30,5])

env = camis.AnisotropicMap(z, hiRes, hexRes,\
                               offset, occupancy_radius, tracking_error)
plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, ax1 = plt.subplots(figsize=(5, 4),nrows = 1, ncols = 1, constrained_layout=True)
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*env.hexSlopeMap, 
                 cmap="nipy_spectral",s=16.0)
cbar = fig.colorbar(cc, ax=ax1,shrink=0.6)
#cc.set_clim(0,50.0)
cbar.set_label('Steepness (deg)')


#cc.set_clim(0,50.0)
#cbar.set_label('Steepness (deg)')





hiRes_elevationMap = np.loadtxt(\
                        open("data/umaRescueArea/UMARescueArea_1mDEM.csv",\
                             "rb"), delimiter=" ", skiprows=0)
offset = np.loadtxt(\
                        open("data/umaRescueArea/UMARescueArea_1mOffset.csv",\
                             "rb"), delimiter=" ", skiprows=0)

#fig = go.Figure(data=[go.Surface(contours = {"z": {"show": True, "size": 0.01, "color":"white"}},\
#                                 z=hiRes_elevationMap[:70,60:140], colorscale = 'haline')])
#
#scene=dict(camera=dict(eye=dict(x=1.15, y=1.15, z=0.8)), #the default values are 1.25, 1.25, 1.25
#           xaxis=dict(),
#           yaxis=dict(),
#           zaxis=dict(),
#           aspectmode='data', #this string can be 'data', 'cube', 'auto', 'manual'
#           #a custom aspectratio is defined as follows:
#           aspectratio=dict(x=1, y=1, z=1)
#           )
#
#fig.update_layout(autosize=True,
#                  width=1024, height=512,
#                  margin=dict(l=0, r=0, b=0, t=0), scene = scene)
#
##fig.show(renderer="iframe")
#fig.write_image("elevation_map.pdf")

hiRes = 1.0
occupancy_radius = 0.5
tracking_error = 0.5

print('TEST_DEMO: DEM is loaded')



#r = 3
#y,x = np.ogrid[-r: r+1, -r: r+1]
#convMatrix = x**2+y**2 <= r**2
#convMatrix = convMatrix.astype(float)
#convMatrix = convMatrix/convMatrix.sum()
#for i in range(5):
#    hiRes_elevationMap = signal.convolve2d(hiRes_elevationMap, \
#                                           convMatrix, \
#                                           mode='same', boundary='symm')

##elevationMap = hiRes_elevationMap[:70,60:140]
##posA = np.asarray([6,10])
##posB = np.asarray([70,10])
##posC = np.asarray([70,50])
##hexRes = 1.0
#
#elevationMap = hiRes_elevationMap[:50,80:140]
#posA = np.asarray([25,5])
#posB = np.asarray([10,40])
#posC = np.asarray([50,20])
#hexRes = 1.0
#    
##elevationMap = hiRes_elevationMap[-70:,-80:]
##posA = np.asarray([45,10])
##posB = np.asarray([70,60])
##posC = np.asarray([10,60])
##hexRes = 1.0
#    
##elevationMap = hiRes_elevationMap[20:60,140:]
##posA = np.asarray([10,5])
##posB = np.asarray([50,30])
##posC = np.asarray([50,5])
##hexRes = 1.0
#
#XX,YY = np.meshgrid(range(elevationMap.shape[1]), range(elevationMap.shape[0]))
##elevationMap = elevationMap - 0.1*YY
##elevationMap = elevationMap - 0.1*XX - 0.1*YY
#
#env = camis.AnisotropicMap(elevationMap, hiRes, hexRes,\
#                               offset, occupancy_radius, tracking_error)

#plt.style.use('default')
#plt.rcParams["font.family"] = "Constantia"
#plt.rcParams['mathtext.fontset'] = 'cm'
#plt.rcParams['mathtext.rm'] = 'serif'
#fig, ax1 = plt.subplots(figsize=(5, 4),nrows = 1, ncols = 1, constrained_layout=True)
#cc = ax1.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*env.hexSlopeMap, 
#                 cmap="nipy_spectral",s=16.0)
#cbar = fig.colorbar(cc, ax=ax1,shrink=0.6)
#cc.set_clim(0,50.0)
#cbar.set_label('Steepness (deg)')

def computeAllPlannings(anisoMapList):
    anisoMapList[0].executeBiPlanning(posB,posA)
#    anisoMapList[1].executeBiPlanning(posC,posB)
#    anisoMapList[2].executeBiPlanning(posA,posC)
#    anisoMapList[3].executeBiPlanning(posC,posA)
#    anisoMapList[4].executeBiPlanning(posB,posC)
    anisoMapList[1].executeBiPlanning(posA,posB)

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



with open("data/sim01/cuadriga_aniso_01.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_01 = camis.CamisDrivingModel(cuadriga_data)
aniso_01.showCAMIS(25)
env_CUAD01_scene01, env_isoCUAD01_scene01 = getMapLists(aniso_01)
computeAllPlannings(env_CUAD01_scene01)
computeAllPlannings(env_isoCUAD01_scene01)

with open("data/sim01/cuadriga_aniso_02.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_02 = camis.CamisDrivingModel(cuadriga_data)
aniso_02.showCAMIS(25)
env_CUAD02_scene01, env_isoCUAD02_scene01 = getMapLists(aniso_02)
computeAllPlannings(env_CUAD02_scene01)
computeAllPlannings(env_isoCUAD02_scene01)

with open("data/sim01/cuadriga_aniso_03.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_03 = camis.CamisDrivingModel(cuadriga_data)
aniso_03.showCAMIS(25)
env_CUAD03_scene01, env_isoCUAD03_scene01 = getMapLists(aniso_03)
computeAllPlannings(env_CUAD03_scene01)
computeAllPlannings(env_isoCUAD03_scene01)

with open("data/sim01/cuadriga_aniso_04.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_04 = camis.CamisDrivingModel(cuadriga_data)
aniso_04.showCAMIS(25)
env_CUAD04_scene01, env_isoCUAD04_scene01 = getMapLists(aniso_04)
computeAllPlannings(env_CUAD04_scene01)
computeAllPlannings(env_isoCUAD04_scene01)

#with open("data/sim01/cuadriga_aniso_05.yml", 'r') as file:
#    cuadriga_data = yaml.full_load(file)
#aniso_05 = camis.CamisDrivingModel(cuadriga_data)
#aniso_05.showCAMIS(55)
#env_CUAD05_scene01, env_isoCUAD05_scene01 = getMapLists(aniso_05)
#computeAllPlannings(env_CUAD05_scene01)
#computeAllPlannings(env_isoCUAD05_scene01)
#
#with open("data/sim01/cuadriga_aniso_06.yml", 'r') as file:
#    cuadriga_data = yaml.full_load(file)
#aniso_06 = camis.CamisDrivingModel(cuadriga_data)
#aniso_06.showCAMIS(55)
#env_CUAD06_scene01, env_isoCUAD06_scene01 = getMapLists(aniso_06)
#computeAllPlannings(env_CUAD06_scene01)
#computeAllPlannings(env_isoCUAD06_scene01)
#
#with open("data/sim01/cuadriga_aniso_07.yml", 'r') as file:
#    cuadriga_data = yaml.full_load(file)
#aniso_07 = camis.CamisDrivingModel(cuadriga_data)
#env_CUAD07_scene01, env_isoCUAD07_scene01 = getMapLists(aniso_07)
#computeAllPlannings(env_CUAD07_scene01)
#computeAllPlannings(env_isoCUAD07_scene01)
#
#with open("data/sim01/cuadriga_aniso_08.yml", 'r') as file:
#    cuadriga_data = yaml.full_load(file)
#aniso_08 = camis.CamisDrivingModel(cuadriga_data)
#env_CUAD08_scene01, env_isoCUAD08_scene01 = getMapLists(aniso_08)
#computeAllPlannings(env_CUAD08_scene01)
#computeAllPlannings(env_isoCUAD08_scene01)
#
#with open("data/sim01/cuadriga_aniso_09.yml", 'r') as file:
#    cuadriga_data = yaml.full_load(file)
#aniso_09 = camis.CamisDrivingModel(cuadriga_data)
#env_CUAD09_scene01, env_isoCUAD09_scene01 = getMapLists(aniso_09)
#computeAllPlannings(env_CUAD09_scene01)
#computeAllPlannings(env_isoCUAD09_scene01)
#
#
#env.show3dDEM()
#
plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, (ax1,ax2) = plt.subplots(figsize=(5, 4),nrows = 1, ncols = 2, constrained_layout=True)
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*env.hexSlopeMap, 
                 cmap="nipy_spectral",s=4.0)
cbar = fig.colorbar(cc, ax=ax1,shrink=0.6)
cbar.set_label('Steepness (deg)')

fig, ax = plt.subplots(constrained_layout=True)
gradient = np.linspace(0,25,26)
slipRatio = 0.07 * np.exp(0.1*gradient)
slipFactorParallel = 1.0 / (1 - slipRatio)
slipAngle = 1.32 * np.exp(0.16*gradient)
slipFactorPerp = 1 / np.cos(slipAngle*3.1416/180.0)
#ax.plot(gradient, slipRatio)
#ax.plot(gradient, slipAngle)
ax.plot(gradient, slipFactorParallel)
ax.plot(gradient, slipFactorPerp)

def showAnisoPath(mapList, color, ax1, ax2, mode):
    for i,anisomap in enumerate(mapList):
        if i < 1:
            anisomap.showPath(fig,ax1,color,mode)
        else:
            anisomap.showPath(fig,ax2,color,mode)
def showIsoPath(mapList, color, ax1, ax2, mode):
    for i,anisomap in enumerate(mapList):
        if i < 1:
            anisomap.showPath(fig,ax1,color,mode)
        else:
            anisomap.showPath(fig,ax2,color,mode)
            
plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, axes = plt.subplots(figsize=(8, 8), \
      nrows = 1, ncols = 3, \
      sharex = 'all', sharey = 'all')
ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]


plt.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.5, top = 1.0, wspace = 0.01, hspace = 0.05)
fig.text(0.5, 0.005, 'X-axis [m]', ha='center')
fig.text(0.005, 0.5, 'Y-axis [m]', va='center', rotation='vertical')

showAnisoPath(env_CUAD01_scene01, 'r', ax1, ax2,'solid')
showIsoPath(env_isoCUAD01_scene01, 'r', ax3, ax3,'dashed')
showAnisoPath(env_CUAD02_scene01, 'c', ax1, ax2,'solid')
showIsoPath(env_isoCUAD02_scene01, 'c', ax3, ax3,'dashed')
showAnisoPath(env_CUAD03_scene01, 'lime', ax1, ax2,'solid')
showIsoPath(env_isoCUAD03_scene01, 'lime', ax3, ax3,'dashed')
showAnisoPath(env_CUAD04_scene01, 'b', ax1, ax2,'solid')
showIsoPath(env_isoCUAD04_scene01, 'b', ax3, ax3,'dashed')
#showAnisoPath(env_CUAD05_scene01, 'm', ax1, ax2,'dashed')
#showIsoPath(env_isoCUAD05_scene01, 'm', ax3, ax3,'dashed')
#showAnisoPath(env_CUAD06_scene01, 'g', ax1, ax2,'dashed')
#showIsoPath(env_isoCUAD06_scene01, 'g', ax3, ax3,'dashed')
#showAnisoPath(env_CUAD07_scene01, 'orange', ax1, ax2,'dotted')
#showIsoPath(env_isoCUAD07_scene01, 'orange', ax3, ax3,'dotted')
#showAnisoPath(env_CUAD08_scene01, 'm', ax1, ax2,'dotted')
#showIsoPath(env_isoCUAD08_scene01, 'm', ax3, ax3,'dotted')   
#showAnisoPath(env_CUAD09_scene01, 'k', ax1, ax2,'dotted')
#showIsoPath(env_isoCUAD09_scene01, 'k', ax3, ax3,'dotted')
#ax1.text(0.5, 0.95, 'Go traverses \n (CAMIS models)', horizontalalignment='center', \
#         verticalalignment='center', transform=ax1.transAxes, fontsize = 12,\
#         color = 'white')
#ax2.text(0.5, 0.95, 'Return traverses \n (CAMIS models)', horizontalalignment='center', \
#         verticalalignment='center', transform=ax2.transAxes, fontsize = 12,\
#         color = 'white')
#ax3.text(0.5, 0.95, 'Go and Return traverses \n (isotropic equivalent models)', horizontalalignment='center', \
#         verticalalignment='center', transform=ax3.transAxes, fontsize = 12,\
#         color = 'white')
##ax1.legend()
#     
for ax in axes:
    cc = ax.scatter(env.hexXmap, env.hexYmap, c = env.hexElevationMap, cmap = cm.gist_earth,s=20)
    ax.scatter(posA[0], posA[1], facecolor = 'r', edgecolor='black', s=60)
    ax.scatter(posB[0], posB[1], facecolor = 'r', edgecolor='black', s=60)
    ax.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
    ax.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
    ax.set_aspect('equal')
fig.tight_layout()



#
#env.executeBiPlanning(posB,posA)
#
#plt.style.use('default')
#plt.rcParams["font.family"] = "Constantia"
#plt.rcParams['mathtext.fontset'] = 'cm'
#plt.rcParams['mathtext.rm'] = 'serif'
#fig, ax = plt.subplots(figsize=(8, 8), \
#      nrows = 1, ncols = 1, \
#      sharex = 'all', sharey = 'all')
#plt.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.5, top = 1.0, wspace = 0.01, hspace = 0.05)
#fig.text(0.5, 0.005, 'X-axis [m]', ha='center')
#fig.text(0.005, 0.5, 'Y-axis [m]', va='center', rotation='vertical')
#
#
#env.showPath(fig,ax,'r','solid')
#
#cc = ax.scatter(env.hexXmap, env.hexYmap, c = env.hexElevationMap, cmap = cm.gist_earth,s=20)
#ax.scatter(posA[0], posA[1], facecolor = 'r', edgecolor='black', s=60)
#ax.scatter(posB[0], posB[1], facecolor = 'r', edgecolor='black', s=60)
#ax.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
#ax.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
#ax.set_aspect('equal')
#fig.tight_layout()


################ ENERGY PLOT ##################
    
        
coeffsLabels = ['CUAD01','CUAD02','CUAD03','CSLI01']
anisoTotalCost = [0,0,0,0]
isoTotalCost = [0,0,0,0]
anisoTR = [0,0]
isoTR = [0,0]
for i in range(2):
    anisoTotalCost[0] = anisoTotalCost[0] + env_CUAD01_scene01[i].pathComputedTotalCost[-1]/3600.0
    anisoTotalCost[1] = anisoTotalCost[1] + env_CUAD02_scene01[i].pathComputedTotalCost[-1]/3600.0
    anisoTotalCost[2] = anisoTotalCost[2] + env_CUAD03_scene01[i].pathComputedTotalCost[-1]/3600.0
    anisoTotalCost[3] = anisoTotalCost[3] + env_CUAD04_scene01[i].pathComputedTotalCost[-1]/3600.0
    isoTotalCost[0] = isoTotalCost[0] + env_isoCUAD01_scene01[i].pathComputedTotalCost[-1]/3600.0
    isoTotalCost[1] = isoTotalCost[1] + env_isoCUAD02_scene01[i].pathComputedTotalCost[-1]/3600.0
    isoTotalCost[2] = isoTotalCost[2] + env_isoCUAD03_scene01[i].pathComputedTotalCost[-1]/3600.0
    isoTotalCost[3] = isoTotalCost[3] + env_isoCUAD04_scene01[i].pathComputedTotalCost[-1]/3600.0
    anisoTR[0] = anisoTR[0] + env_CUAD02_scene01[i].pathComputedTotalCostwithRisk[-1]/3600.0
    anisoTR[1] = anisoTR[1] + env_CUAD03_scene01[i].pathComputedTotalCostwithRisk[-1]/3600.0
    isoTR[0] = isoTR[0] + env_isoCUAD02_scene01[i].pathComputedTotalCostwithRisk[-1]/3600.0
    isoTR[1] = isoTR[1] + env_isoCUAD03_scene01[i].pathComputedTotalCostwithRisk[-1]/3600.0

x = np.arange(len(coeffsLabels))  # the label locations
x2 = np.arange(2)+1
width = 0.4  # the width of the bars

fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
rects3 = ax.bar(x2 - 0.45/2, anisoTR, 0.45, label='Isotropic (ρ = 0.8)', color='lime')
rects4 = ax.bar(x2 + 0.45/2, isoTR, 0.45, label='Isotropic (ρ = 0.8)', color='g')
rects1 = ax.bar(x - width/2, anisoTotalCost, width, label='Isotropic (ρ = 0.8)', color='r')
rects2 = ax.bar(x + width/2, isoTotalCost, width, label='Isotropic (ρ = 0.8)', color = 'b')

for i,rect in enumerate(rects1):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T) + '\nAh',
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, -3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='top',color = 'w')
for i,rect in enumerate(rects2):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T) + '\nAh',
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, -3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='top',color = 'w')
for i,rect in enumerate(rects2):
        isoT = rect.get_height()
        anisoT = rects1[i].get_height()
        gain = (isoT - anisoT)/isoT * 100
        ax.annotate('Gain = ' + '{0:.2f}'.format(gain) + '%',
                    xy=(rect.get_x(), isoT),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
for i,rect in enumerate(rects3):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T) + '\nAh',
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom',color = 'k')
for i,rect in enumerate(rects4):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T) + '\nAh',
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom',color = 'k')
for i,rect in enumerate(rects4):
        isoT = rect.get_height()
        anisoT = rects3[i].get_height()
        gain = (isoT - anisoT)/isoT * 100
        ax.annotate('Gain = ' + '{0:.2f}'.format(gain) + '%',
                    xy=(rect.get_x(), isoT),
                    xytext=(0, 25),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')   
        
ax.grid(True, which='both')   
#autolabel(rects1)
#autolabel(rects2)
ax.set_ylabel('Total Cost [Ah]')
ax.set_xlabel('CAMIS')
ax.set_xticks(x)
ax.set_xticklabels(coeffsLabels)
ax.legend(('Anisotropic Total Cost + Risk','Isotropic Total Cost + Risk', 'Anisotropic Total Cost','Isotropic Total Cost'))
plt.minorticks_on()  
plt.show()















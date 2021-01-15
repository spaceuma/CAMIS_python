# -*- coding: utf-8 -*-
"""
Numerical Test for CAMIS paper
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from context import camis
import copy
import zipfile
import yaml
from matplotlib.ticker import FormatStrFormatter
try:
    from scipy import signal
except:
    raise ImportError('ERROR: scipy module could not be imported')
import plotly.graph_objects as go
from matplotlib.colors import LightSource

# =============================================================================
## LOADING DEM
# =============================================================================

#v = x.*exp(-x.^2-y.^2-z.^2);

xx = np.linspace(-3.0,3.0,50)
yy = np.linspace(-4.5,4.5,50)
xmesh, ymesh = np.meshgrid(xx,yy)

z = np.ones_like(xmesh)
for j,y in enumerate(yy):
    for i,x in enumerate(xx):
        z[j,i] =  3*(1-x)**2*np.exp(-(x**2)/1.0 - (y+1)**2/1.0) - \
        10*(x/5 - x**3 - y**5)*np.exp(-x**2-y**2) - \
        1/3*np.exp(-(x+1)**2 - y**2) 

xmesh *= 10.0
ymesh *= 10.0
z = z/9.0

plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, ax1 = plt.subplots(figsize=(3, 4),nrows = 1, ncols = 1, constrained_layout=True)
cc = ax1.scatter(xmesh, ymesh, c = z, 
                 cmap=cm.gist_earth,s=16.0)
cbar = fig.colorbar(cc, ax=ax1,shrink=0.6, location = 'top')
cbar.set_label('Elevation [deg]')
ax1.set_xlabel('X [m]')
ax1.set_ylabel('Y [m]')
ax1.set_aspect('equal')

hiRes = 1.0
hexRes = 0.5
offset = (0,0)
occupancy_radius = 0.5
tracking_error = 0.5
#posA = np.asarray([20,45])
#posB = np.asarray([30,5])
posA = np.asarray([27,5])
posB = np.asarray([20,45])
env = camis.AnisotropicMap(z, hiRes, hexRes,\
                               offset, occupancy_radius, tracking_error)
plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, ax1 = plt.subplots(figsize=(3.2, 3.9),nrows = 1, ncols = 1, constrained_layout=True)
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*env.hexSlopeMap, 
                 cmap="nipy_spectral",s=16.0, vmin = 0.0, vmax = 25.0)
cbar = fig.colorbar(cc, ax=ax1,shrink=0.6, location = 'top', \
                    ticks=[0,5,10,15,20,25])
#cc.set_clim(0,50.0)
cbar.set_label('Steepness α [deg]')
ax1.set_xlim([0,48.0])
ax1.set_ylim([0,48.0])
ax1.set_xlabel('X [m]')
ax1.set_ylabel('Y [m]')
ax1.set_aspect('equal')


plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, ax1 = plt.subplots(figsize=(3.2, 3.9),nrows = 1, ncols = 1, constrained_layout=True)
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*np.arctan2(env.hexAspectMap[1],env.hexAspectMap[0]),
                 cmap="gist_rainbow",s=16.0, vmin = -180.0, vmax = 180.0)
cbar = fig.colorbar(cc, ax=ax1,shrink=0.6, location = 'top', ticks=[-180, -90, 0, 90, 180])
cbar.set_label('Aspect Angle [deg]')
ax1.set_xlim([0,48.0])
ax1.set_ylim([0,48.0])
ax1.set_xlabel('X [m]')
ax1.set_ylabel('Y [m]')
ax1.set_aspect('equal')


#cc.set_clim(0,50.0)
#cbar.set_label('Steepness (deg)')





hiRes_elevationMap = np.loadtxt(\
                        open("data/umaRescueArea/UMARescueArea_1mDEM.csv",\
                             "rb"), delimiter=" ", skiprows=0)
offset = np.loadtxt(\
                        open("data/umaRescueArea/UMARescueArea_1mOffset.csv",\
                             "rb"), delimiter=" ", skiprows=0)



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

with open("data/sim01/cuadriga_aniso_05.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_05 = camis.CamisDrivingModel(cuadriga_data)
aniso_05.showCAMIS(25)
env_CUAD05_scene01, env_isoCUAD05_scene01 = getMapLists(aniso_05)
computeAllPlannings(env_CUAD05_scene01)
computeAllPlannings(env_isoCUAD05_scene01)

with open("data/sim01/cuadriga_aniso_06.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_06 = camis.CamisDrivingModel(cuadriga_data)
aniso_06.showCAMIS(25)
env_CUAD06_scene01, env_isoCUAD06_scene01 = getMapLists(aniso_06)
computeAllPlannings(env_CUAD06_scene01)
computeAllPlannings(env_isoCUAD06_scene01)

plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, ax = plt.subplots(figsize=(5,3),constrained_layout=True)
p1 = aniso_01.showModelData('anisotropy',fig,ax,'r','solid',25)
p2 = aniso_02.showModelData('anisotropy',fig,ax,'m','solid',25)
p3 = aniso_03.showModelData('anisotropy',fig,ax,'orange','solid',25)
p4 = aniso_04.showModelData('anisotropy',fig,ax,'b','solid',25)
p5 = aniso_05.showModelData('anisotropy',fig,ax,'c','solid',25)
p6 = aniso_06.showModelData('anisotropy',fig,ax,'g','solid',25)
p7, = ax.plot([0], marker='None',
           linestyle='None', label='Model A')
ax.set_xlabel('Steepness α [degrees]')
ax.set_ylabel('Anisotropy')
l1 = ax.legend([p7,p1,p2,p3,p7,p7,p4,p5,p6], [r'$Model A$'] + \
               ['$ρ = 0.3$', '$ρ = 0.6$', \
                '$ρ = 0.9$'] + ['']  + \
                [r'$Model B$'] + \
                [r'$ρ = 0.3$', r'$ρ = 0.6$', \
                r'$ρ = 0.9$'], ncol = 2)
plt.show()


def plotModels(model, fig, axes,label,modelname):
    model.showModelData('ascent-cost',fig,axes,'m','dashed',25)
    model.showModelData('lateral-cost',fig,axes,'g','dashed',25)
    model.showModelData('descent-cost',fig,axes,'b','dashed',25)
    axes.legend(('$C(α,\pm π)$','$C(α,\pm π/2)$','$C(α, 0)$'), loc = 2)
#    axes.set_ylabel(label)
    axes.text(0.05, 0.35, modelname, horizontalalignment='left', \
         verticalalignment='bottom', transform=axes.transAxes, fontsize = 10,\
         color = 'k')
plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig,(axes1, axes2, axes3) = plt.subplots(figsize=(6, 6), \
      nrows = 3, ncols = 2)
plotModels(aniso_01,fig,axes1[0],'Power/Speed [As/m]','$ρ = 0.3$')
plotModels(aniso_02,fig,axes2[0],'Power/Speed [As/m]','$ρ = 0.6$')
plotModels(aniso_03,fig,axes3[0],'Power/Speed [As/m]','$ρ = 0.9$')
plotModels(aniso_04,fig,axes1[1],'Power/Speed [As/m]','$ρ = 0.3$')
plotModels(aniso_05,fig,axes2[1],'Power/Speed [As/m]','$ρ = 0.6$')
plotModels(aniso_06,fig,axes3[1],'Power/Speed [As/m]','$ρ = 0.9$')
for ax in axes1:
#    ax.set_ylim([0,300])
    ax.grid(True, which='both') 
for ax in axes2:
#    ax.set_ylim([0,3300])
    ax.grid(True, which='both')
for ax in axes3:
#    ax.set_ylim([0,110])
    ax.grid(True, which='both')
axes1[0].set_title('Model A')
axes1[1].set_title('Model B')
plt.subplots_adjust(left = 0.1, right = 0.99, bottom = 0.075, top = 0.95, \
                    wspace = 0.15, hspace = 0.15)
fig.text(0.5, 0.0075, 'Steepness α [deg]', ha='center')
axes2[0].set_ylabel('Power/Speed [As/m]')
plt.minorticks_on() 
#axes.set_xlabel('Steepness [deg]')

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
showAnisoPath(env_CUAD02_scene01, 'm', ax1, ax2,'solid')
showIsoPath(env_isoCUAD02_scene01, 'm', ax3, ax3,'dashed')
showAnisoPath(env_CUAD03_scene01, 'orange', ax1, ax2,'solid')
showIsoPath(env_isoCUAD03_scene01, 'orange', ax3, ax3,'dashed')
showAnisoPath(env_CUAD04_scene01, 'b', ax1, ax2,'solid')
showIsoPath(env_isoCUAD04_scene01, 'b', ax3, ax3,'dashed')
showAnisoPath(env_CUAD05_scene01, 'c', ax1, ax2,'dashed')
showIsoPath(env_isoCUAD05_scene01, 'c', ax3, ax3,'dashed')
showAnisoPath(env_CUAD06_scene01, 'g', ax1, ax2,'dashed')
showIsoPath(env_isoCUAD06_scene01, 'g', ax3, ax3,'dashed')
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
plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
coeffsLabels = ['Model A\n$ρ= 0.3$','Model A\n$ρ = 0.6$',\
                'Model A\n$ρ = 0.9$','Model B\n$ρ = 0.3$',\
                'Model B\n$ρ = 0.6$','Model B\n$ρ = 0.9$']
anisoTotalCost = [0,0,0,0,0,0]
isoTotalCost = [0,0,0,0,0,0]
anisoTR = [0,0]
isoTR = [0,0]
for i in range(2):
    anisoTotalCost[0] = anisoTotalCost[0] + env_CUAD01_scene01[i].pathComputedTotalCost[-1]#/3600.0
    anisoTotalCost[1] = anisoTotalCost[1] + env_CUAD02_scene01[i].pathComputedTotalCost[-1]#/3600.0
    anisoTotalCost[2] = anisoTotalCost[2] + env_CUAD03_scene01[i].pathComputedTotalCost[-1]#/3600.0
    anisoTotalCost[3] = anisoTotalCost[3] + env_CUAD04_scene01[i].pathComputedTotalCost[-1]#/3600.0
    anisoTotalCost[4] = anisoTotalCost[4] + env_CUAD05_scene01[i].pathComputedTotalCost[-1]#/3600.0
    anisoTotalCost[5] = anisoTotalCost[5] + env_CUAD06_scene01[i].pathComputedTotalCost[-1]#/3600.0
    isoTotalCost[0] = isoTotalCost[0] + env_isoCUAD01_scene01[i].pathComputedTotalCost[-1]#/3600.0
    isoTotalCost[1] = isoTotalCost[1] + env_isoCUAD02_scene01[i].pathComputedTotalCost[-1]#/3600.0
    isoTotalCost[2] = isoTotalCost[2] + env_isoCUAD03_scene01[i].pathComputedTotalCost[-1]#/3600.0
    isoTotalCost[3] = isoTotalCost[3] + env_isoCUAD04_scene01[i].pathComputedTotalCost[-1]#/3600.0
    isoTotalCost[4] = isoTotalCost[4] + env_isoCUAD05_scene01[i].pathComputedTotalCost[-1]#/3600.0
    isoTotalCost[5] = isoTotalCost[5] + env_isoCUAD06_scene01[i].pathComputedTotalCost[-1]#/3600.0
#    anisoTR[0] = anisoTR[0] + env_CUAD02_scene01[i].pathComputedTotalCostwithRisk[-1]/3600.0
#    anisoTR[1] = anisoTR[1] + env_CUAD03_scene01[i].pathComputedTotalCostwithRisk[-1]/3600.0
#    isoTR[0] = isoTR[0] + env_isoCUAD02_scene01[i].pathComputedTotalCostwithRisk[-1]/3600.0
#    isoTR[1] = isoTR[1] + env_isoCUAD03_scene01[i].pathComputedTotalCostwithRisk[-1]/3600.0

x = np.arange(len(coeffsLabels))  # the label locations
x2 = np.arange(2)+1
width = 0.4  # the width of the bars

fig, ax = plt.subplots(figsize=(7.5,4.5), constrained_layout=True)
#rects3 = ax.bar(x2 - 0.45/2, anisoTR, 0.45, label='Isotropic (ρ = 0.8)', color='lime')
#rects4 = ax.bar(x2 + 0.45/2, isoTR, 0.45, label='Isotropic (ρ = 0.8)', color='g')
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
#for i,rect in enumerate(rects3):
#        T = rect.get_height()
#        ax.annotate('{0:.2f}'.format(T) + '\nAh',
#                    xy=(rect.get_x() + rect.get_width() / 2, T),
#                    xytext=(0, 3),  # 3 points vertical offset
#                    textcoords="offset points",
#                    ha='center', va='bottom',color = 'k')
#for i,rect in enumerate(rects4):
#        T = rect.get_height()
#        ax.annotate('{0:.2f}'.format(T) + '\nAh',
#                    xy=(rect.get_x() + rect.get_width() / 2, T),
#                    xytext=(0, 3),  # 3 points vertical offset
#                    textcoords="offset points",
#                    ha='center', va='bottom',color = 'k')
#for i,rect in enumerate(rects4):
#        isoT = rect.get_height()
#        anisoT = rects3[i].get_height()
#        gain = (isoT - anisoT)/isoT * 100
#        ax.annotate('Gain = ' + '{0:.2f}'.format(gain) + '%',
#                    xy=(rect.get_x(), isoT),
#                    xytext=(0, 25),  # 3 points vertical offset
#                    textcoords="offset points",
#                    ha='center', va='bottom')   
        
ax.grid(True, which='both')   
#autolabel(rects1)
#autolabel(rects2)
ax.set_ylabel('Total Cost [As]')
#ax.set_xlabel('CAMIS')
ax.set_xticks(x)
ax.set_xticklabels(coeffsLabels)
ax.legend(('Anisotropic Total Cost','Isotropic Total Cost'))
plt.minorticks_on()  
plt.show()


####3d PATHS###
plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig = plt.figure(figsize=(13.2, 4.5))
ax1 = fig.add_subplot(1, 3, 1, projection='3d')
ax2 = fig.add_subplot(1, 3, 2, projection='3d')
ax3 = fig.add_subplot(1, 3, 3, projection='3d')

XX = np.zeros_like(z)
YY = np.zeros_like(z)
for i in range(50):
    for j in range(50):
        XX[j,i] = i
        YY[j,i] = j
ls = LightSource(270, 45)
rgb = ls.shade(z, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
surf1 = ax1.plot_surface(XX,YY,z, rstride=1, cstride=1, facecolors=rgb,\
                       linewidth=0, antialiased=False, shade=False)
surf2 = ax2.plot_surface(XX,YY,z, rstride=1, cstride=1, facecolors=rgb,\
                       linewidth=0, antialiased=False, shade=False)
surf3 = ax3.plot_surface(XX,YY,z, rstride=1, cstride=1, facecolors=rgb,\
                       linewidth=0, antialiased=False, shade=False)
ax1.plot(env_CUAD01_scene01[0].path[:,0], env_CUAD01_scene01[0].path[:,1], \
        env_CUAD01_scene01[0].pathElevation,linestyle='dashed',color = 'r')
ax1.plot(env_CUAD02_scene01[0].path[:,0], env_CUAD02_scene01[0].path[:,1], \
        env_CUAD02_scene01[0].pathElevation,linestyle='dashed',color = 'm')
ax1.plot(env_CUAD03_scene01[0].path[:,0], env_CUAD03_scene01[0].path[:,1], \
        env_CUAD03_scene01[0].pathElevation,linestyle='dashed',color = 'orange')
ax1.plot(env_CUAD04_scene01[0].path[:,0], env_CUAD04_scene01[0].path[:,1], \
        env_CUAD04_scene01[0].pathElevation,linestyle='dashed',color = 'b')
ax1.plot(env_CUAD05_scene01[0].path[:,0], env_CUAD05_scene01[0].path[:,1], \
        env_CUAD05_scene01[0].pathElevation,linestyle='dashed',color = 'c')
ax1.plot(env_CUAD06_scene01[0].path[:,0], env_CUAD06_scene01[0].path[:,1], \
        env_CUAD06_scene01[0].pathElevation,linestyle='dashed',color = 'lime')
ax1.set_xlabel('X [m]')
ax1.set_ylabel('Y [m]')
ax1.set_zlabel('Z [m]')
ax1.set_zticks([-0.6, 0.1, 0.8])
ax1.view_init(78.0, -150.0)
ax1.set_facecolor('w')
ax2.plot(env_CUAD01_scene01[1].path[:,0], env_CUAD01_scene01[1].path[:,1], \
        env_CUAD01_scene01[1].pathElevation,linestyle='dashed',color = 'r')
ax2.plot(env_CUAD02_scene01[1].path[:,0], env_CUAD02_scene01[1].path[:,1], \
        env_CUAD02_scene01[1].pathElevation,linestyle='dashed',color = 'm')
ax2.plot(env_CUAD03_scene01[1].path[:,0], env_CUAD03_scene01[1].path[:,1], \
        env_CUAD03_scene01[1].pathElevation,linestyle='dashed',color = 'orange')
ax2.plot(env_CUAD04_scene01[1].path[:,0], env_CUAD04_scene01[1].path[:,1], \
        env_CUAD04_scene01[1].pathElevation,linestyle='dashed',color = 'b')
ax2.plot(env_CUAD05_scene01[1].path[:,0], env_CUAD05_scene01[1].path[:,1], \
        env_CUAD05_scene01[1].pathElevation,linestyle='dashed',color = 'c')
ax2.plot(env_CUAD06_scene01[1].path[:,0], env_CUAD06_scene01[1].path[:,1], \
        env_CUAD06_scene01[1].pathElevation,linestyle='dashed',color = 'lime')
ax2.set_xlabel('X [m]')
ax2.set_ylabel('Y [m]')
ax2.set_zlabel('Z [m]')
ax2.set_zticks([-0.6, 0.1, 0.8])
ax2.view_init(78.0, -150.0)
ax2.set_facecolor('w')

ax3.plot(env_isoCUAD01_scene01[0].path[:,0], env_isoCUAD01_scene01[0].path[:,1], \
         env_isoCUAD01_scene01[0].pathElevation,linestyle='dashed',color = 'r')
ax3.plot(env_isoCUAD02_scene01[0].path[:,0], env_isoCUAD02_scene01[0].path[:,1], \
         env_isoCUAD02_scene01[0].pathElevation,linestyle='dashed',color = 'm')
ax3.plot(env_isoCUAD03_scene01[0].path[:,0], env_isoCUAD03_scene01[0].path[:,1], \
         env_isoCUAD03_scene01[0].pathElevation,linestyle='dashed',color = 'orange')
ax3.plot(env_isoCUAD04_scene01[0].path[:,0], env_isoCUAD04_scene01[0].path[:,1], \
         env_isoCUAD04_scene01[0].pathElevation,linestyle='dashed',color = 'b')
ax3.plot(env_isoCUAD05_scene01[0].path[:,0], env_isoCUAD05_scene01[0].path[:,1], \
         env_isoCUAD05_scene01[0].pathElevation,linestyle='dashed',color = 'c')
ax3.plot(env_isoCUAD06_scene01[0].path[:,0], env_isoCUAD06_scene01[0].path[:,1], \
         env_isoCUAD06_scene01[0].pathElevation,linestyle='dashed',color = 'lime')
ax3.set_xlabel('X [m]')
ax3.set_ylabel('Y [m]')
ax3.set_zlabel('Z [m]')
ax3.set_zticks([-0.6, 0.1, 0.8])
ax3.view_init(78.0, -150.0)
ax3.set_facecolor('w')
fig.tight_layout()

#MAX = 6
#for direction in (-1, 1):
#    for point in np.diag(direction * MAX * np.array([1,1,1])):
#        ax.plot([point[0]+12.5], [point[1]+12.5], [point[2]+2.0], 'w')
        
##### COMPUTATIONAL TIMES ###
anisoTime = [0,0,0,0]
isoTime = [0,0,0,0]
for i in range(2):
    print('Aniso')
    print(env_CUAD01_scene01[i].elapsedTime)
    print(env_CUAD02_scene01[i].elapsedTime)
    print(env_CUAD03_scene01[i].elapsedTime)
    print(env_CUAD04_scene01[i].elapsedTime)
    print(env_CUAD05_scene01[i].elapsedTime)
    print(env_CUAD06_scene01[i].elapsedTime)
    print('Iso')
    print(env_isoCUAD01_scene01[i].elapsedTime)
    print(env_isoCUAD02_scene01[i].elapsedTime)
    print(env_isoCUAD03_scene01[i].elapsedTime)
    print(env_isoCUAD04_scene01[i].elapsedTime)
    print(env_isoCUAD05_scene01[i].elapsedTime)
    print(env_isoCUAD06_scene01[i].elapsedTime)
    anisoTime[0] = anisoTime[0] + env_CUAD01_scene01[i].elapsedTime
    anisoTime[1] = anisoTime[1] + env_CUAD02_scene01[i].elapsedTime
    anisoTime[2] = anisoTime[2] + env_CUAD03_scene01[i].elapsedTime
    anisoTime[3] = anisoTime[3] + env_CUAD04_scene01[i].elapsedTime
    isoTime[0] = isoTime[0] + env_isoCUAD01_scene01[i].elapsedTime
    isoTime[1] = isoTime[1] + env_isoCUAD02_scene01[i].elapsedTime
    isoTime[2] = isoTime[2] + env_isoCUAD03_scene01[i].elapsedTime 
    isoTime[3] = isoTime[3] + env_isoCUAD04_scene01[i].elapsedTime     

print(anisoTime)
print(isoTime)


fig = go.Figure(data=[go.Surface(contours = {"z": {"show": True, "size": 0.01, "color":"white"}},\
                                 z=z, colorscale = 'haline'), \
    go.Scatter3d(
    x=env_CUAD01_scene01[0].path[:,0], y=env_CUAD01_scene01[0].path[:,1], \
    z=env_CUAD01_scene01[0].pathElevation,
    marker=dict(
        size=2,
        color='red'
    )), \
    go.Scatter3d(
    x=env_CUAD01_scene01[1].path[:,0], y=env_CUAD01_scene01[1].path[:,1], \
    z=env_CUAD01_scene01[1].pathElevation,
    marker=dict(
        size=2,
        color='red'
    )), \
    go.Scatter3d(
    x=env_CUAD02_scene01[0].path[:,0], y=env_CUAD02_scene01[0].path[:,1], \
    z=env_CUAD02_scene01[0].pathElevation,
    marker=dict(
        size=2,
        color='magenta'
    )), \
    go.Scatter3d(
    x=env_CUAD02_scene01[1].path[:,0], y=env_CUAD02_scene01[1].path[:,1], \
    z=env_CUAD02_scene01[1].pathElevation,
    marker=dict(
        size=2,
        color='magenta'
    )), \
    go.Scatter3d(
    x=env_CUAD03_scene01[0].path[:,0], y=env_CUAD03_scene01[0].path[:,1], \
    z=env_CUAD03_scene01[0].pathElevation,
    marker=dict(
        size=2,
        color='orange'
    )), \
    go.Scatter3d(
    x=env_CUAD03_scene01[1].path[:,0], y=env_CUAD03_scene01[1].path[:,1], \
    z=env_CUAD03_scene01[1].pathElevation,
    marker=dict(
        size=2,
        color='orange'
    )), \
    go.Scatter3d(
    x=env_CUAD04_scene01[0].path[:,0], y=env_CUAD04_scene01[0].path[:,1], \
    z=env_CUAD04_scene01[0].pathElevation,
    marker=dict(
        size=2,
        color='blue'
    )), \
    go.Scatter3d(
    x=env_CUAD04_scene01[1].path[:,0], y=env_CUAD04_scene01[1].path[:,1], \
    z=env_CUAD04_scene01[1].pathElevation,
    marker=dict(
        size=2,
        color='blue'
    )), \
    go.Scatter3d(
    x=env_CUAD05_scene01[0].path[:,0], y=env_CUAD05_scene01[0].path[:,1], \
    z=env_CUAD05_scene01[0].pathElevation,
    marker=dict(
        size=2,
        color='cyan'
    )), \
    go.Scatter3d(
    x=env_CUAD05_scene01[1].path[:,0], y=env_CUAD05_scene01[1].path[:,1], \
    z=env_CUAD05_scene01[1].pathElevation,
    marker=dict(
        size=2,
        color='cyan'
    )), \
    go.Scatter3d(
    x=env_CUAD06_scene01[0].path[:,0], y=env_CUAD06_scene01[0].path[:,1], \
    z=env_CUAD06_scene01[0].pathElevation,
    marker=dict(
        size=2,
        color='green'
    )), \
    go.Scatter3d(
    x=env_CUAD06_scene01[1].path[:,0], y=env_CUAD06_scene01[1].path[:,1], \
    z=env_CUAD06_scene01[1].pathElevation,
    marker=dict(
        size=2,
        color='green'
    )) 
])

scene=dict(camera=dict(up=dict(x=0, y=0, z=1),\
                       center=dict(x=0, y=0, z=-0.2),\
                       eye=dict(x=1.5, y=-1.5, z=0.7)), #the default values are 1.25, 1.25, 1.25
           xaxis=dict(),
           yaxis=dict(),
           zaxis=dict(),
           aspectmode='data', #this string can be 'data', 'cube', 'auto', 'manual'
           #a custom aspectratio is defined as follows:
           aspectratio=dict(x=1, y=1, z=1)
           )

fig.update_layout(autosize=False,
                  width=2048, height=512,
                  margin=dict(l=0, r=0, b=0, t=0), scene = scene)
#fig.show(renderer="iframe")
fig.write_image("elevation_map.pdf")















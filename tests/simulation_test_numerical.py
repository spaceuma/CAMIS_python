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
                 cmap="nipy_spectral",s=16.0, vmin = 0.0, vmax = 25.0, rasterized=True)
cbar = fig.colorbar(cc, ax=ax1,shrink=0.6, location = 'top', \
                    ticks=[0,5,10,15,20,25])
#cc.set_clim(0,50.0)
cbar.set_label('Steepness α [deg]')
ax1.set_xlim([0,48.0])
ax1.set_ylim([0,48.0])
ax1.set_xlabel('X [m]')
ax1.set_ylabel('Y [m]')
ax1.set_aspect('equal')
plt.savefig('numerical_tests_steepnessMap_reduced.pdf',dpi=300)


plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, ax1 = plt.subplots(figsize=(3.2, 3.9),nrows = 1, ncols = 1, constrained_layout=True)
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*np.arctan2(env.hexAspectMap[1],env.hexAspectMap[0]),
                 cmap="gist_rainbow",s=16.0, vmin = -180.0, vmax = 180.0, rasterized=True)
cbar = fig.colorbar(cc, ax=ax1,shrink=0.6, location = 'top', ticks=[-180, -90, 0, 90, 180])
cbar.set_label('Aspect Angle [deg]')
ax1.set_xlim([0,48.0])
ax1.set_ylim([0,48.0])
ax1.set_xlabel('X [m]')
ax1.set_ylabel('Y [m]')
ax1.set_aspect('equal')
plt.savefig('numerical_tests_aspectMap_reduced.pdf',dpi=300)


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
fig, ax = plt.subplots(figsize=(5,2.5),constrained_layout=True)
p1 = aniso_01.showModelData('anisotropy',fig,ax,'r','solid',25)
p2 = aniso_02.showModelData('anisotropy',fig,ax,'m','solid',25)
p3 = aniso_03.showModelData('anisotropy',fig,ax,'orange','solid',25)
p4 = aniso_04.showModelData('anisotropy',fig,ax,'b','solid',25)
p5 = aniso_05.showModelData('anisotropy',fig,ax,'c','solid',25)
p6 = aniso_06.showModelData('anisotropy',fig,ax,'lime','solid',25)
p7, = ax.plot([0], marker='None',
           linestyle='None', label='Wheel Model')
ax.set_xlabel('Steepness α [degrees]', fontsize = 12)
ax.set_ylabel('Anisotropy', fontsize = 12)
ax.text(-1.8,5.7,'ϒ',{'family':'DejaVu Sans'},rotation = 'vertical', fontsize = 12)
l1 = ax.legend([p7,p1,p2,p3,p7,p7,p4,p5,p6], [r'$Wheel \ Model$'] + \
               ['$ρ = 0.3$', '$ρ = 0.6$', \
                '$ρ = 0.9$'] + ['']  + \
                [r'$Track \ Model$'] + \
                [r'$ρ = 0.3$', r'$ρ = 0.6$', \
                r'$ρ = 0.9$'], ncol = 2, fontsize = 10)
ax.set_xlim([0,25])
plt.minorticks_on()
plt.grid(b=True,which='minor', linestyle = '--')
plt.grid(b=True,which='major', linewidth = 1)
plt.show()


def plotModels(model, fig, axes,label,modelname):
    axes.plot(0,0,alpha = 0.0)
    model.showModelData('ascent-cost',fig,axes,'m','dashed',25)
    model.showModelData('lateral-cost',fig,axes,'g','dashed',25)
    model.showModelData('descent-cost',fig,axes,'b','dashed',25)
    model.showModelData('nominal-cost',fig,axes,'orange','dashed',25)
    axes.legend((modelname,'$C(α, β = \pm π)$','$C(α, β = \pm π/2)$',\
                 '$C(α, β = 0)$', '$C_n(α)$'), loc = 2, ncol = 2, fontsize = 10)
#    axes.set_ylabel(label)
#    axes.text(0.05, 0.35, modelname, horizontalalignment='left', \
#         verticalalignment='bottom', transform=axes.transAxes, fontsize = 10,\
#         color = 'k')
plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig,(axes1, axes2, axes3) = plt.subplots(figsize=(7, 6), \
      nrows = 3, ncols = 2)
plotModels(aniso_01,fig,axes1[0],'Power/Speed [As/m]','$\mathbf{ρ = 0.3}$')
axes1[0].grid(b=True,which='minor', linestyle = '--')
axes1[0].grid(b=True,which='major', linewidth = 1)
axes1[0].minorticks_on() 
plotModels(aniso_02,fig,axes2[0],'Power/Speed [As/m]','$\mathbf{ρ = 0.6}$')
axes2[0].grid(b=True,which='minor', linestyle = '--')
axes2[0].grid(b=True,which='major', linewidth = 1)
axes2[0].minorticks_on() 
plotModels(aniso_03,fig,axes3[0],'Power/Speed [As/m]','$\mathbf{ρ = 0.9}$')
axes3[0].grid(b=True,which='minor', linestyle = '--')
axes3[0].grid(b=True,which='major', linewidth = 1)
axes3[0].minorticks_on() 
plotModels(aniso_04,fig,axes1[1],'Power/Speed [As/m]','$\mathbf{ρ = 0.3}$')
axes1[1].grid(b=True,which='minor', linestyle = '--')
axes1[1].grid(b=True,which='major', linewidth = 1)
axes1[1].minorticks_on() 
plotModels(aniso_05,fig,axes2[1],'Power/Speed [As/m]','$\mathbf{ρ = 0.6}$')
axes2[1].grid(b=True,which='minor', linestyle = '--')
axes2[1].grid(b=True,which='major', linewidth = 1)
axes2[1].minorticks_on() 
plotModels(aniso_06,fig,axes3[1],'Power/Speed [As/m]','$\mathbf{ρ = 0.9}$')
axes3[1].grid(b=True,which='minor', linestyle = '--')
axes3[1].grid(b=True,which='major', linewidth = 1)
axes3[1].minorticks_on() 
for ax in (axes1, axes2, axes3):
#    ax.set_ylim([0,300])
    ax[0].grid(True, which='both')
    ax[1].grid(True, which='both')
    ax[0].set_xlim([0,25])
    ax[1].set_xlim([0,25])
#    ax[0].set_ylim([0,9.0])
    ax[1].set_ylim([0, 25.0])
axes1[0].set_title('Wheel Model', fontsize = 14)
axes1[1].set_title('Track Model', fontsize = 14)
plt.subplots_adjust(left = 0.07, right = 0.99, bottom = 0.075, top = 0.95, \
                    wspace = 0.15, hspace = 0.15)
fig.text(0.5, 0.0075, 'Steepness α [deg]', ha='center', fontsize = 14)
axes2[0].set_ylabel('Energy per Distance [J/m]', fontsize = 14)
plt.grid(b=True,which='minor', linestyle = '--')
plt.grid(b=True,which='major', linewidth = 1)
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
coeffsLabels = ['Wheel\nModel\n$ρ= 0.3$','Wheel\nModel\n$ρ = 0.6$',\
                'Wheel\nModel\n$ρ = 0.9$','Track\nModel\n$ρ = 0.3$',\
                'Track\nModel\n$ρ = 0.6$','Track\nModel\n$ρ = 0.9$']
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
width = 0.48  # the width of the bars

fig, ax = plt.subplots(figsize=(7,3), constrained_layout=True)
#rects3 = ax.bar(x2 - 0.45/2, anisoTR, 0.45, label='Isotropic (ρ = 0.8)', color='lime')
#rects4 = ax.bar(x2 + 0.45/2, isoTR, 0.45, label='Isotropic (ρ = 0.8)', color='g')
rects1 = ax.bar(x - width/2, anisoTotalCost, width, label='Isotropic (ρ = 0.8)', color='r')
rects2 = ax.bar(x + width/2, isoTotalCost, width, label='Isotropic (ρ = 0.8)', color = 'b')

for i,rect in enumerate(rects1):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T),
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, -3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='top',color = 'w', fontsize = 12)
for i,rect in enumerate(rects2):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T),
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, -3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='top',color = 'w', fontsize = 12)
for i,rect in enumerate(rects2):
        isoT = rect.get_height()
        anisoT = rects1[i].get_height()
        gain = (isoT - anisoT)/isoT * 100
        ax.annotate('{0:.2f}'.format(gain) + '%',
                    xy=(rect.get_x(), isoT),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize = 12)
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
ax.set_ylabel('Total Cost [J]', fontsize = 12)
ax.set_ylim([0,1000])
#ax.set_xlabel('CAMIS')
ax.set_xticks(x)
ax.set_xticklabels(coeffsLabels, fontsize = 12)
ax.legend(('Anisotropic','Isotropic'), fontsize = 12)
plt.minorticks_on()  
plt.show()


####3d PATHS###
XX = np.zeros_like(z)
YY = np.zeros_like(z)
for i in range(50):
    for j in range(50):
        XX[j,i] = i
        YY[j,i] = j
        
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = env.hexElevationMap, \
                cmap = cm.gist_earth,s=20,vmin = -1.0, vmax = 1.0)

plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig = plt.figure(figsize=(13.2, 4.5))
ax1 = fig.add_subplot(1, 3, 1, projection='3d')
ax2 = fig.add_subplot(1, 3, 2, projection='3d')
ax3 = fig.add_subplot(1, 3, 3, projection='3d')


ls = LightSource(270, 45)
rgb = ls.shade(z, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
surf1 = ax1.plot_surface(XX,YY,z, rstride=1, cstride=1, facecolors=rgb,\
                       linewidth=0, antialiased=False, shade=False,\
                       vmin = -1.0, vmax = 1.0, rasterized=True)
surf2 = ax2.plot_surface(XX,YY,z, rstride=1, cstride=1, facecolors=rgb,\
                       linewidth=0, antialiased=False, shade=False,\
                       vmin = -1.0, vmax = 1.0, rasterized=True)
surf3 = ax3.plot_surface(XX,YY,z, rstride=1, cstride=1, facecolors=rgb,\
                       linewidth=0, antialiased=False, shade=False,\
                       vmin = -1.0, vmax = 1.0, rasterized=True)
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
ax1.set_xlabel('X [m]',labelpad=-8)
ax1.set_ylabel('Y [m]',labelpad=-8)
ax1.set_zlabel('Z [m]',labelpad=-6)
ax1.tick_params(axis="x",direction="in", pad=-5)
ax1.tick_params(axis="y",direction="in", pad=-5)
ax1.tick_params(axis="z",direction="in", pad=-2)
ax1.set_zticks([-0.6, 0.1, 0.8])
ax1.view_init(78.0, -150.0)
ax1.set_xlim([0,50])
ax1.set_ylim([0,50])
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
ax2.set_xlabel('X [m]',labelpad=-8)
ax2.set_ylabel('Y [m]',labelpad=-8)
ax2.set_zlabel('Z [m]',labelpad=-6)
ax2.tick_params(axis="x",direction="in", pad=-5)
ax2.tick_params(axis="y",direction="in", pad=-5)
ax2.tick_params(axis="z",direction="in", pad=-2)
ax2.set_zticks([-0.6, 0.1, 0.8])
ax2.view_init(78.0, -150.0)
ax2.set_xlim([0,50])
ax2.set_ylim([0,50])
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
ax3.set_xlabel('X [m]',labelpad=-8)
ax3.set_ylabel('Y [m]',labelpad=-8)
ax3.set_zlabel('Z [m]',labelpad=-6)
ax3.tick_params(axis="x",direction="in", pad=-5)
ax3.tick_params(axis="y",direction="in", pad=-5)
ax3.tick_params(axis="z",direction="in", pad=-2)
ax3.set_zticks([-0.6, 0.1, 0.8])
ax3.view_init(78.0, -150.0)
ax3.set_xlim([0,50])
ax3.set_ylim([0,50])
ax3.set_facecolor('w')

ax1.text(posA[0],posA[1],0.2,'$X_o$',color='w', size=15)
ax2.text(posA[0],posA[1],0.2,'$X_o$',color='w', size=15)
ax3.text(posA[0],posA[1],0.2,'$X_o$',color='w', size=15)
ax1.text(posB[0],posB[1],0.2,'$X_g$',color='w', size=15)
ax2.text(posB[0],posB[1],0.2,'$X_g$',color='w', size=15)
ax3.text(posB[0],posB[1],0.2,'$X_g$',color='w', size=15)

ax1.text(50,35,0.8,'$Anisotropic$\n $X_o \Rightarrow X_g$',color='k', size=14)
ax2.text(50,35,0.8,'$Anisotropic$\n $X_g \Rightarrow X_o$',color='k', size=14)
ax3.text(50,35,0.8,'$Isotropic$\n $Both \ ways$',color='k', size=14)


fig.tight_layout()
plt.subplots_adjust(left = 0.01, right = 1.0, bottom = 0.0, top = 1.0, \
                    wspace = 0.0, hspace = 0.0)
cbar_ax = fig.add_axes([0.01, 0.05, 0.15, 0.025])
legend_ax1 = fig.add_axes([0.1, 0.18, 0.3, 0.1])
legend_ax2 = fig.add_axes([0.1, 0.18, 0.63, 0.1])
legend_ax1.plot(0,0, alpha = 0, label = '$\mathbf{Wheel \ Model}$')
legend_ax1.plot(0,0,'r', linestyle = 'dashed', label='$ρ = 0.3$')
legend_ax1.plot(0,0,'m', linestyle = 'dashed', label='$ρ = 0.6$')
legend_ax1.plot(0,0,'orange', linestyle = 'dashed', label='$ρ = 0.9$')
legend_ax1.grid(b=False)
legend_ax1.axis('off')
legend_ax1.legend(fontsize=12)
legend_ax2.plot(0,0, alpha = 0, label = '$\mathbf{Track \ Model}$')
legend_ax2.plot(0,0,'b', linestyle = 'dashed', label='$ρ = 0.3$')
legend_ax2.plot(0,0,'c', linestyle = 'dashed', label='$ρ = 0.6$')
legend_ax2.plot(0,0,'lime', linestyle = 'dashed', label='$ρ = 0.9$')
legend_ax2.grid(b=False)
legend_ax2.axis('off')
legend_ax2.legend(fontsize=12)
#plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
#plt.grid(b=False)
#plt.axis('off')
cbar = fig.colorbar(cc, cax=cbar_ax, orientation = 'horizontal')
cbar.set_alpha(1.0)
cbar.set_label('Elevation [m]',color='k', labelpad=-40,x = 0.4,fontsize=14)
cbar.ax.tick_params(colors='k', size = 2)
cbar.outline.set_visible(True)
#cbar.outline.set_edgecolor('k')
cbar.ax.yaxis.set_tick_params(color='k')
cbar.set_ticks([-1.0, -0.5, 0.0, 0.5, 1.0])
#cbar.outline.set_linewidth(0.5)
plt.savefig('numerical_tests_paths_reduced.pdf',dpi=300)
#MAX = 6
#for direction in (-1, 1):
#    for point in np.diag(direction * MAX * np.array([1,1,1])):
#        ax1.plot([point[0]+12.5], [point[1]+12.5], [point[2]+2.0], 'w')
#plt.show()

        
##### COMPUTATIONAL TIMES ###
coeffsLabels = ['Wheel\nModel\n$ρ= 0.3$','Wheel\nModel\n$ρ = 0.6$',\
                'Wheel\nModel\n$ρ = 0.9$','Track\nModel\n$ρ = 0.3$',\
                'Track\nModel\n$ρ = 0.6$','Track\nModel\n$ρ = 0.9$']
anisoTime = [0,0,0,0,0,0]
isoTime = [0,0,0,0,0,0]
for i in range(2):
    anisoTime[0] = anisoTime[0] + env_CUAD01_scene01[i].elapsedTime
    anisoTime[1] = anisoTime[1] + env_CUAD02_scene01[i].elapsedTime
    anisoTime[2] = anisoTime[2] + env_CUAD03_scene01[i].elapsedTime
    anisoTime[3] = anisoTime[3] + env_CUAD04_scene01[i].elapsedTime
    anisoTime[4] = anisoTime[4] + env_CUAD05_scene01[i].elapsedTime
    anisoTime[5] = anisoTime[5] + env_CUAD06_scene01[i].elapsedTime
    isoTime[0] = isoTime[0] + env_isoCUAD01_scene01[i].elapsedTime
    isoTime[1] = isoTime[1] + env_isoCUAD02_scene01[i].elapsedTime
    isoTime[2] = isoTime[2] + env_isoCUAD03_scene01[i].elapsedTime 
    isoTime[3] = isoTime[3] + env_isoCUAD04_scene01[i].elapsedTime   
    isoTime[4] = isoTime[4] + env_isoCUAD05_scene01[i].elapsedTime
    isoTime[5] = isoTime[5] + env_isoCUAD06_scene01[i].elapsedTime

x = np.arange(len(coeffsLabels))  # the label locations
x2 = np.arange(2)+1
width = 0.4  # the width of the bars

fig, ax = plt.subplots(figsize=(7,3), constrained_layout=True)
#rects3 = ax.bar(x2 - 0.45/2, anisoTR, 0.45, label='Isotropic (ρ = 0.8)', color='lime')
#rects4 = ax.bar(x2 + 0.45/2, isoTR, 0.45, label='Isotropic (ρ = 0.8)', color='g')
rects1 = ax.bar(x - width/2, anisoTime, width, label='Isotropic (ρ = 0.8)', color='r')
rects2 = ax.bar(x + width/2, isoTime, width, label='Isotropic (ρ = 0.8)', color = 'b')

for i,rect in enumerate(rects1):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T),
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, -3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='top',color = 'w', fontsize = 12)
for i,rect in enumerate(rects2):
        T = rect.get_height()
        ax.annotate('{0:.2f}'.format(T),
                    xy=(rect.get_x() + rect.get_width() / 2, T),
                    xytext=(0, -3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='top',color = 'w', fontsize = 12)
for i,rect in enumerate(rects1):
        anisoT = rect.get_height()
        isoT = rects2[i].get_height()
        gain = (isoT - anisoT)/isoT * 100
        ax.annotate('{0:.2f}'.format(gain) + '%',
                    xy=(rects2[i].get_x(), anisoT),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize = 12)
ax.grid(True, which='both')   
#autolabel(rects1)
#autolabel(rects2)
ax.set_ylabel('Computational Time [s]', fontsize = 12)
#ax.set_xlabel('CAMIS')
ax.set_xticks(x)
ax.set_xticklabels(coeffsLabels, fontsize = 12)
ax.set_ylim([0,40])
ax.legend(('Anisotropic','Isotropic'), fontsize = 12)
plt.minorticks_on()  
plt.show()


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















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


#cc.set_clim(0,50.0)
#cbar.set_label('Steepness (deg)')



hiRes = 1.0
occupancy_radius = 0.5
tracking_error = 0.5

print('TEST_DEMO: DEM is loaded')

def computeAllPlannings(anisoMapList):
    anisoMapList[0].executeBiPlanning(posB,posA)
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



with open("data/sim01/cuadriga_aniso_03.yml", 'r') as file:
    cuadriga_data = yaml.full_load(file)
aniso_01 = camis.CamisDrivingModel(cuadriga_data)
#aniso_01.showCAMIS(25)
env_CUAD01_scene01, env_isoCUAD01_scene01 = getMapLists(aniso_01)

env_isoCUAD01_scene01[0].executeBiPlanning(posB,posA)



# computeAllPlannings(env_CUAD01_scene01)
# computeAllPlannings(env_isoCUAD01_scene01)

# def showAnisoPath(mapList, color, ax1, ax2, mode):
#     for i,anisomap in enumerate(mapList):
#         if i < 1:
#             anisomap.showPath(fig,ax1,color,mode)
#         else:
#             anisomap.showPath(fig,ax2,color,mode)
# def showIsoPath(mapList, color, ax1, ax2, mode):
#     for i,anisomap in enumerate(mapList):
#         if i < 1:
#             anisomap.showPath(fig,ax1,color,mode)
#         else:
#             anisomap.showPath(fig,ax2,color,mode)
            
# plt.style.use('default')
# plt.rcParams["font.family"] = "Constantia"
# plt.rcParams['mathtext.fontset'] = 'cm'
# plt.rcParams['mathtext.rm'] = 'serif'
# fig, axes = plt.subplots(figsize=(8, 8), \
#       nrows = 1, ncols = 3, \
#       sharex = 'all', sharey = 'all')
# ax1 = axes[0]
# ax2 = axes[1]
# ax3 = axes[2]


# plt.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.5, top = 1.0, wspace = 0.01, hspace = 0.05)
# fig.text(0.5, 0.005, 'X-axis [m]', ha='center')
# fig.text(0.005, 0.5, 'Y-axis [m]', va='center', rotation='vertical')

# showAnisoPath(env_CUAD01_scene01, 'r', ax1, ax2,'solid')
# showIsoPath(env_isoCUAD01_scene01, 'r', ax3, ax3,'dashed')

# for ax in axes:
#     cc = ax.scatter(env.hexXmap, env.hexYmap, c = env.hexElevationMap, cmap = cm.gist_earth,s=20)
#     ax.scatter(posA[0], posA[1], facecolor = 'r', edgecolor='black', s=60)
#     ax.scatter(posB[0], posB[1], facecolor = 'r', edgecolor='black', s=60)
#     ax.set_xlim([env.xMap[0,2], env.xMap[-1,-4]])
#     ax.set_ylim([env.yMap[0,0], env.yMap[-1,-1]])
#     ax.set_aspect('equal')
# fig.tight_layout()
# #
# env_isoCUAD01_scene01[0].showResults()
# env_CUAD01_scene01[0].showResults()
# #
# #env_isoCUAD01_scene01[1].showResults()
# #env_CUAD01_scene01[1].showResults()
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 10:37:29 2020

@author: Richi
"""

######################## SHOW PATHS ##################
def showAnisoPath(mapList, color, ax1, ax2, mode):
    for i,anisomap in enumerate(mapList):
        if i < 4:
            anisomap.showPath(fig,ax1,color,mode)
        else:
            anisomap.showPath(fig,ax2,color,mode)
def showIsoPath(mapList, color, ax1, ax2, mode):
    for i,anisomap in enumerate(mapList):
        if i < 4:
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
showIsoPath(env_isoCUAD01_scene01, 'r', ax3, ax3,'solid')
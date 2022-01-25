# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 16:34:53 2022

@author: Richi
"""

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
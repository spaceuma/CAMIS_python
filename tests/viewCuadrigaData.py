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


# Load previous_execution/202002processedcuadrigadata on spyder

import numpy as np

import data.cuadrigaData.cuadriga_reader as cr
from context import camis
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yaml
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import scipy.signal
from scipy.optimize import curve_fit

def fittingSlip(x,x1):
    return np.piecewise(x, 
                        [x < 0, x >= 0],
                        [lambda x: - x1 * x, lambda x: x1 * x])



#    return np.piecewise(x, 
#                        [x < 0, x >= 0],
#                        [lambda x: x2 * x + Co, lambda x: x1 * x + Co])
 
bounds = ([0.0],[1.0])
slipArray = 1.0 - 1.0/speedCommand*np.asarray(speedList)
poptSlip,_ = curve_fit(fittingSlip, (signedGradient), slipArray, bounds = bounds)




def fittingCost(x,Kmg,rho):
    return np.piecewise(x, 
                        [x < 0, x >= 0],
                        [lambda x: Kmg*(rho - np.tan(np.pi/180*x))/(1 - poptSlip[0]*x),\
                         lambda x: Kmg*(rho - np.tan(np.pi/180*x))/(1 - poptSlip[0]*x)])
          
bounds = ([0.0,0.0],[np.inf,1.0])
popt,_ = curve_fit(fittingCost, (signedGradient), costList, bounds = bounds)

plt.style.use('seaborn-darkgrid')
fig2 = plt.figure(figsize=(7, 4))
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
left, width = 0.08, 0.75
bottom, height = 0.12, 0.65
spacing = 0.01
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.15, height]
ax2 = plt.axes(rect_scatter)
ax2.set_facecolor('xkcd:light blue grey')
#ax2.set_title('$\longleftarrow$ Ascent Direction      Descent Direction $\longrightarrow$',y = -0.13)
plt.grid('True')
ax2.scatter(signedGradient, costList, color='r', s=1, rasterized=True)
ax2.set_ylim(0,np.max(costList))
ax2.set_xlim(-signedGradient.max(),signedGradient.max())
ax2.plot(np.arange(-0.0, 16.0), fittingCost(np.arange(-0.0, 16.0), popt[0], popt[1]),'b', linestyle = 'solid')
ax2.plot(np.arange(-16.0, 1.0), fittingCost(np.arange(-16.0, 1.0), popt[0], popt[1]),'b', linestyle = 'dashed')
ax2.legend(['$0.45 - tan_θ$', '$0.45 + tan_θ$','Experimental Data'], fontsize=12, facecolor = 'xkcd:pale grey', frameon = True)
plt.xlabel('Pitch θ [degrees]', fontsize = 14)
plt.ylabel('Energy per Distance [As/m]', fontsize = 14)
ax_histx = plt.axes(rect_histx)
ax_histx.tick_params(direction='in', labelbottom=False)
#ax_histx.set_facecolor('xkcd:light grey blue')
ax_histy = plt.axes(rect_histy)
ax_histy.tick_params(direction='in', labelleft=False)
#ax_histy.set_facecolor('xkcd:light grey blue')
    
binwidth = 0.25
binsX = np.arange(-signedGradient.max(), signedGradient.max() + binwidth, binwidth)
binsY = np.arange(0, np.max(costList) + binwidth, binwidth)
ax_histx.hist(signedGradient, bins=binsX, color = "m", rasterized=True)
ax_histy.hist(costList, bins=binsY, orientation='horizontal', color = "m", rasterized=True)
    
ax_histx.set_xlim(ax2.get_xlim())
ax_histy.set_ylim(ax2.get_ylim())
fig2.tight_layout()
plt.savefig('ascentDescentCuadrigaCostFunction_reduced.pdf',dpi=300)




fig2 = plt.figure(figsize=(7, 4))
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
left, width = 0.08, 0.75
bottom, height = 0.12, 0.65
spacing = 0.01
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.15, height]
ax2 = plt.axes(rect_scatter)
ax2.set_facecolor('xkcd:light grey blue')
#ax2.set_title('$\longleftarrow$ Ascent Direction      Descent Direction $\longrightarrow$',y = -0.13)
plt.grid('True')
ax2.scatter(signedGradient, slipArray, color='r', s=1, rasterized=True)
ax2.set_ylim(-0.5,0.5)
ax2.set_xlim(-signedGradient.max(),signedGradient.max())
ax2.plot(np.arange(-16.0, 16.0), fittingSlip(np.arange(-16.0, 16.0), poptSlip[0]),'g', linestyle = 'solid')
#ax2.plot(np.arange(-16.0, 1.0), fittingSlip(np.arange(-16.0, 1.0), poptSlip[0], poptSlip[1]),'g', linestyle = 'dashed')
ax2.legend(['$s_r = 3^{-19}$','Experimental Data'], fontsize=12, facecolor = 'xkcd:pale grey', frameon = True)
plt.xlabel('Pitch θ [degrees]', fontsize = 14)
plt.ylabel('Slip Ratio', fontsize = 14)
ax2.yaxis.set_label_coords(-0.05,0.5)
ax_histx = plt.axes(rect_histx)
ax_histx.tick_params(direction='in', labelbottom=False)
ax_histy = plt.axes(rect_histy)
ax_histy.tick_params(direction='in', labelleft=False)
    
binwidth = 0.25
binsX = np.arange(-signedGradient.max(), signedGradient.max() + binwidth, binwidth)
binsY = np.arange(np.min(slipArray), np.max(slipArray) + 0.0025, 0.0025)
ax_histx.hist(signedGradient, bins=binsX, color = "m", rasterized=True)
ax_histy.hist(slipArray, bins=binsY, orientation='horizontal', color = "m", rasterized=True)
    
ax_histx.set_xlim(ax2.get_xlim())
ax_histy.set_ylim(ax2.get_ylim())
fig2.tight_layout()
plt.savefig('ascentDescentCuadrigaSlipFunction_reduced.pdf',dpi=300)
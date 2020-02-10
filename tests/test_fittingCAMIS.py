# -*- coding: utf-8 -*-
"""
DEMO script: building CAMIS based on Cuadriga experimental data
@author: J.Ricardo Sanchez Ibanez (ricardosan@uma.es)

This script is a reference demo to understand how CAMIS is built using data
obtained by a mobile robot.

"""

import numpy as np

import data.cuadrigaData.cuadriga_reader as cr
from context import camis
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yaml
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import scipy.signal
from scipy.optimize import curve_fit

fittingMode = 'DEM'
#fittingMode = 'IMU'


#if fittingMode == 'DEM':
    
hiRes_elevationMap = np.loadtxt(open("data/terrainData/UMATerrain_10cmDEM.csv",\
                                         "rb"), delimiter=" ", skiprows=0)
#hiRes_posX = np.loadtxt(open("data/terrainData/UMARescueArea_10cmPosX.csv",\
#                                         "rb"), delimiter=" ", skiprows=0)
#hiRes_posY = np.loadtxt(open("data/terrainData/UMARescueArea_10cmPosY.csv",\
#                                         "rb"), delimiter=" ", skiprows=0)
    
hiRes = 0.1

#offset = np.loadtxt(\
#                        open("data/terrainData/UMATerrain_10cmOffset.csv",\
#                             "rb"), delimiter=" ", skiprows=0)
#
#with open("data/cuadriga.yml", 'r') as file:
#    cuadriga_data = yaml.full_load(file)
#r1 = camis.CamisDrivingModel(cuadriga_data)
#
#env = camis.AnisotropicMap(hiRes_elevationMap[0:2600,0:1500], hiRes, 0.4,\
#                               offset)
#env.computeVecCostMap(r1)


#env = camis.PDEM(hiRes_elevationMap, hiRes, .5, (hiRes_posX[0,0],hiRes_posY[0,0]))
#env.smoothMap(1.0) #Tracking error is not considered here



"01 - The CSV files collected during preliminar experiments are loaded"
csvFiles = ['data/cuadrigaData/20190129/2019_01_29_10_53_36.txt',
            'data/cuadrigaData/20190129/2019_01_29_10_59_06.txt',
            'data/cuadrigaData/20190129/2019_01_29_11_05_32.txt',
            'data/cuadrigaData/20190129/2019_01_29_11_08_35.txt',
            'data/cuadrigaData/20190129/2019_01_29_11_16_00.txt',
            'data/cuadrigaData/20190129/2019_01_29_11_20_23.txt',
            'data/cuadrigaData/20190129/2019_01_29_11_24_09.txt',
#            'data/cuadrigaData/20190531/JornadasRescate01.txt',
#            'data/cuadrigaData/20190531/JornadasRescate02.txt',
#            'data/cuadrigaData/20190531/JornadasRescate03.txt',
            'data/cuadrigaData/20190624/2019_06_24_11_17_55.txt',
            'data/cuadrigaData/20190624/2019_06_24_11_22_21.txt',
            'data/cuadrigaData/20190624/2019_06_24_11_32_17.txt',
            'data/cuadrigaData/20190624/2019_06_24_11_36_06.txt',
            'data/cuadrigaData/20190624/2019_06_24_11_52_47.txt',
            'data/cuadrigaData/20190624/2019_06_24_11_57_35.txt'
            ]

#cr.showData('cuadrigaData/20190624/2019_06_24_11_17_55.txt')

"02 - CAMIS input data is extracted"
betaX = []
betaY = []
beta = []
cost = []
gradient = []
speedList = []
imuGradientList = []
betaList = []
sdList = []
costList = []

fig1, axes1 = plt.subplots(constrained_layout=True)
env.showMap('slope-deg',fig1,axes1)
#axes1.set_xlim(10, 55)
#axes1.set_ylim(60, 110)
#axes1.tick_params(axis='both', which='major', labelsize=12)

fig2 = plt.figure()
ax1 = fig2.add_subplot(311)
ax1.set_ylabel('Cost (As/m)', fontsize = 14)
plt.grid(True)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
ax2 = fig2.add_subplot(312)
ax2.set_ylabel('Roll (degrees)', fontsize = 14)
plt.grid(True)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
ax3 = fig2.add_subplot(313)
ax3.set_ylabel('Pitch (degrees)', fontsize = 14)
ax3.set_xlabel('Traversed distance (m)', fontsize = 14)
plt.grid(True)
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
       
ax1.get_shared_x_axes().join(ax1, ax2)
ax1.set_xticklabels([])
ax2.get_shared_x_axes().join(ax2, ax3)
ax2.set_xticklabels([])

ax1.tick_params(axis='both', which='major', labelsize=12)
ax2.tick_params(axis='both', which='major', labelsize=12)
ax3.tick_params(axis='both', which='major', labelsize=12)

for index,file in enumerate(csvFiles):
    " Cost C according to robot orientation Roll-Pitch-Yaw "
    posX, posY, heading, Roll, Pitch, Yaw, Current, Speed, Distance, Segment, GPSspeed, Time = cr.readCuadrigaData(file)
    " Conversion from Roll-Pitch-Yaw to Gradient-Beta_angle "
    if fittingMode == 'IMU':
        G, Bx, By = camis.rpy2ab(Roll,Pitch,Yaw)
        betaX = betaX + Bx
        betaY = betaY + By
        B = np.arctan2(By,Bx)
        beta = beta + B.tolist()
        currentCost = np.divide(Current,Speed)
        cost = cost + currentCost.tolist()
        gradient = gradient + G
        
    if fittingMode == 'DEM':
        imuG, Bx, By = camis.rpy2ab(Roll,Pitch,Yaw)
        G, B = env.getBeta(posX - offset[0], posY - offset[1], heading)
        sd = env.getSD(posX - offset[0], posY - offset[1])
        betaX = np.cos(B)
        betaY = np.sin(B)
        B = np.arctan2(betaY,betaX)
        beta = beta + B.tolist()
        currentCost = np.divide(Current,Speed)
        cost = cost + currentCost.tolist()
        G = np.abs(G)*180/np.pi
        gradient = gradient + G.tolist()
    axes1.plot(posX - offset[0], posY - offset[1])
    
    ax3.plot(Distance,np.convolve(Current, np.ones((10,))/10, mode='same'))
#    Speed = np.convolve(Speed, np.ones((10,))/10, mode='same')
    Speed = scipy.ndimage.median_filter(Speed,size=10)
    Speed = Speed[5:-5]
    ax2.plot(Distance[5:-5],Speed)
    speedList = speedList + Speed.tolist()
    pitch = camis.ab2pitch(G,B)
    imuG = np.convolve(pitch, np.ones((10,))/10, mode='same')
    imuG = imuG[5:-5]
    imuGradientList = imuGradientList + imuG.tolist()
    B = B[5:-5]
    betaList = betaList + B.tolist()
    sdList = sdList + sd
    currentCost = currentCost[5:-5]
    costList = costList + currentCost.tolist()
    ax1.plot(Distance[5:-5],imuG)

fig2.tight_layout()

signedGradient = np.zeros_like(betaList)

for i,b in enumerate(betaList):
    if betaList[i] > np.pi/2 or betaList[i] < -np.pi/2:
        signedGradient[i] = -imuGradientList[i]
    else:
        signedGradient[i] = imuGradientList[i]
        
speedCommand = 0.5

def fittingSlip(x,x1,x2):
    return np.piecewise(x, 
                        [x < 0, x >= 0],
                        [lambda x: - x2 * x, lambda x: x1 * x])



#    return np.piecewise(x, 
#                        [x < 0, x >= 0],
#                        [lambda x: x2 * x + Co, lambda x: x1 * x + Co])
 
bounds = ([0.0,0.0],[1.0,1.0])
slipArray = 1.0 - 1.0/speedCommand*np.asarray(speedList)
poptSlip,_ = curve_fit(fittingSlip, (signedGradient), slipArray, bounds = bounds)

def fittingCost(x,Kmg,rho):
    return np.piecewise(x, 
                        [x < 0, x >= 0],
                        [lambda x: Kmg*(rho - np.tan(np.pi/180*x))/(1 - poptSlip[1]*x),\
                         lambda x: Kmg*(rho - np.tan(np.pi/180*x))/(1 - poptSlip[0]*x)])
          
bounds = ([0.0,0.0],[np.inf,1.0])
popt,_ = curve_fit(fittingCost, (signedGradient), costList, bounds = bounds)


fig2 = plt.figure(figsize=(4, 4))
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
left, width = 0.12, 0.65
bottom, height = 0.12, 0.65
spacing = 0.01
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]
ax2 = plt.axes(rect_scatter)
ax2.set_facecolor('xkcd:pale grey')
#ax2.set_title('$\longleftarrow$ Ascent Direction      Descent Direction $\longrightarrow$',y = -0.13)
plt.grid('True')
ax2.scatter(signedGradient, costList, color='r', s=1)
ax2.set_ylim(0,np.max(costList))
ax2.set_xlim(-signedGradient.max(),signedGradient.max())
ax2.plot(np.arange(-0.0, 16.0), fittingCost(np.arange(-0.0, 16.0), popt[0], popt[1]),'b')
ax2.plot(np.arange(-16.0, 1.0), fittingCost(np.arange(-16.0, 1.0), popt[0], popt[1]),'g')
ax2.legend(['$\mathbb{C}_d$', '$\mathbb{C}_a$','Cost Samples'], fontsize='small')
plt.xlabel('Pitch θ [degrees]')
plt.ylabel('Motors Current Consumption [A]')
ax_histx = plt.axes(rect_histx)
ax_histx.tick_params(direction='in', labelbottom=False)
ax_histy = plt.axes(rect_histy)
ax_histy.tick_params(direction='in', labelleft=False)
    
binwidth = 0.25
binsX = np.arange(-signedGradient.max(), signedGradient.max() + binwidth, binwidth)
binsY = np.arange(0, np.max(costList) + binwidth, binwidth)
ax_histx.hist(signedGradient, bins=binsX, color = "m")
ax_histy.hist(costList, bins=binsY, orientation='horizontal', color = "m")
    
ax_histx.set_xlim(ax2.get_xlim())
ax_histy.set_ylim(ax2.get_ylim())
fig2.tight_layout()


fig2 = plt.figure(figsize=(4, 4))
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
left, width = 0.15, 0.62
bottom, height = 0.12, 0.65
spacing = 0.01
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]
ax2 = plt.axes(rect_scatter)
ax2.set_facecolor('xkcd:pale grey')
#ax2.set_title('$\longleftarrow$ Ascent Direction      Descent Direction $\longrightarrow$',y = -0.13)
plt.grid('True')
ax2.scatter(signedGradient, slipArray, color='r', s=1)
ax2.set_ylim(np.min(slipArray),np.max(slipArray))
ax2.set_xlim(-signedGradient.max(),signedGradient.max())
ax2.plot(np.arange(-0.0, 16.0), fittingSlip(np.arange(-0.0, 16.0), poptSlip[0], poptSlip[1]),'b')
ax2.plot(np.arange(-16.0, 1.0), fittingSlip(np.arange(-16.0, 1.0), poptSlip[0], poptSlip[1]),'g')
ax2.legend(['$s_d |α|$', '$s_a |α|$','Slip Ratio Samples'], fontsize='small')
plt.xlabel('Pitch θ [degrees]')
plt.ylabel('Slip Ratio')
ax_histx = plt.axes(rect_histx)
ax_histx.tick_params(direction='in', labelbottom=False)
ax_histy = plt.axes(rect_histy)
ax_histy.tick_params(direction='in', labelleft=False)
    
binwidth = 0.25
binsX = np.arange(-signedGradient.max(), signedGradient.max() + binwidth, binwidth)
binsY = np.arange(np.min(slipArray), np.max(slipArray) + 0.0025, 0.0025)
ax_histx.hist(signedGradient, bins=binsX, color = "m")
ax_histy.hist(slipArray, bins=binsY, orientation='horizontal', color = "m")
    
ax_histx.set_xlim(ax2.get_xlim())
ax_histy.set_ylim(ax2.get_ylim())
fig2.tight_layout()





fig4, axes4 = plt.subplots()
axes4.scatter(signedGradient, speedList)
axes4.plot(np.arange(-15.0, 16.0), fittingSlip(np.arange(-15.0, 16.0), poptSlip[0], poptSlip[1]),'r')

#fig4 = plt.figure()
#axes4 = fig4.add_subplot(111,projection='3d')
#axes4.scatter(gradient, sdList, cost)


if fittingMode == 'DEM':
    maxCost = 40.0
    gradient_threshold = np.max(gradient)
if fittingMode == 'IMU':
    maxCost = 40.0
    gradient_threshold = np.max(gradient)
   

#beta.append(np.pi)
#cost.append(maxCost)
#gradient.append(gradient_threshold)
#beta.append(0)
#cost.append(maxCost)
#gradient.append(gradient_threshold)

#beta.append(np.pi/2)
#cost.append(maxCost*0.75)
#gradient.append(gradient_threshold/2)
#beta.append(-np.pi/2)
#cost.append(maxCost*0.75)
#gradient.append(gradient_threshold/2)
#beta.append(np.pi/2)
#cost.append(40.0)
#gradient.append(gradient_threshold)
#beta.append(-np.pi/2)
#cost.append(40.0)
#gradient.append(gradient_threshold)

sigma = np.ones_like(gradient)
#sigma[[-2,-1]] = 0.01



"03 - The functions Cd, Ca, Cl1 and Cl2 are obtained through a fitting process"
beta = [0 if np.isnan(x) else x for x in beta]

#CdRoots, CaRoots, Cl1Roots, Cl2Roots = \
#camis.computeDirCosts(gradient,beta,cost,sigma)

r1.fitCAMIS(gradient,beta,cost,sigma)

#drivingRoots = camis.computeDrivingDirCosts(gradient,beta,cost,sigma)


#r1 = camis.CamisModel.fromRoots(CdRoots,CaRoots,Cl1Roots,Cl2Roots, gradient_threshold)
#
#isoRoots = (np.asarray(CaRoots) + np.asarray(CdRoots) + np.asarray(Cl1Roots) + np.asarray(Cl2Roots))/4
#
#rIsoMed = camis.CamisModel.fromRoots(isoRoots,isoRoots,isoRoots,isoRoots, gradient_threshold)
#rIsoAsc = camis.CamisModel.fromRoots(CaRoots,CaRoots,CaRoots,CaRoots, gradient_threshold)

#AniCoLUT = camis.computeAniCoLUT(CdRoots, CaRoots, Cl1Roots, Cl2Roots, gradient_threshold)

"04 - Representation of CAMIS together with raw data"
#camis.showCAMIS(CdRoots, CaRoots, Cl1Roots, Cl2Roots,beta,gradient,cost)
r1.showCAMIS()
#rIsoMed.showCAMIS()
#rIsoAsc.showCAMIS()

"05 - Saving CAMIS"
#if fittingMode == 'DEM':
#    camis.saveCamis('data/camisRoots/cuadriga_camis_dem.csv',CdRoots,CaRoots,Cl1Roots,Cl2Roots,\
#                r1.anicoLUT)
#    camis.saveCamis('data/camisRoots/cuadriga_camis_dem_iso_med.csv',isoRoots,isoRoots,isoRoots,isoRoots,\
#                rIsoMed.anicoLUT)
#    camis.saveCamis('data/camisRoots/cuadriga_camis_dem_iso_asc.csv',CaRoots,CaRoots,CaRoots,CaRoots,\
#                rIsoAsc.anicoLUT)
#if fittingMode == 'IMU':
#    camis.saveCamis('data/camisRoots/cuadriga_camis_imu.csv',CdRoots,CaRoots,Cl1Roots,Cl2Roots,\
#                r1.anicoLUT)
#    camis.saveCamis('data/camisRoots/cuadriga_camis_imu_iso_med.csv',isoRoots,isoRoots,isoRoots,isoRoots,\
#                rIsoMed.anicoLUT)
#    camis.saveCamis('data/camisRoots/cuadriga_camis_imu_iso_asc.csv',CaRoots,CaRoots,CaRoots,CaRoots,\
#                rIsoAsc.anicoLUT)
    

#rDriving = camis.CamisDrivingModel(25.0)
#
#rDriving.fitCAMIS(gradient,beta,cost,sigma)
#
#rDriving.showCAMIS()




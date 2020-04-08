#============================CAMIS library=====================================
#           Continuous Anisotropic Model for Inclined Surfaces
#          Author: J. Ricardo Sanchez Ibanez (ricardosan@uma.es)
# -----------------------------------------------------------------------------
#                           cost_function.py
#   This file contains a library of python functions dedicated to the 
#==============================================================================

import numpy as np
import scipy
from scipy import optimize
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from bisect import bisect_left
import csv
import math

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

deg2rad = np.pi/180
rad2deg = 180/np.pi


# =============================================================================
#    rpy2ab -> Mapping from Tait-Bryan angles to Slope-Traverse angles
# =============================================================================

def rpy2ab(Roll,Pitch,Yaw):
    BetaX = []
    BetaY = []
    aspectVectorX = []
    aspectVectorY = []
    headingVectorX = []
    headingVectorY=[]
    
    # Gradient is computed
    Gradient = []
    for i, r in enumerate(Roll):
        Gradient.append(rad2deg*np.arccos(np.cos(deg2rad*Roll[i])*\
                                          np.cos(deg2rad*Pitch[i])))
        
    for i, r in enumerate(Roll):
        aspect = [np.cos(deg2rad*Yaw[i])*np.sin(deg2rad*Pitch[i])*\
                  np.cos(deg2rad*Roll[i]) + np.sin(deg2rad*Yaw[i])*\
                  np.sin(deg2rad*Roll[i]), np.sin(deg2rad*Yaw[i])*\
                  np.sin(deg2rad*Pitch[i])*np.cos(deg2rad*Roll[i])\
                  - np.cos(deg2rad*Yaw[i])*np.sin(deg2rad*Roll[i])]
        aspect = aspect/np.linalg.norm(aspect)
        aspectVectorX.append(aspect[0])
        aspectVectorY.append(aspect[1])
        heading = [np.cos(deg2rad*Yaw[i])*np.cos(deg2rad*Pitch[i]), \
                   np.sin(deg2rad*Yaw[i])*np.cos(deg2rad*Pitch[i])]
        heading = heading/np.linalg.norm(heading)
        headingVectorX.append(heading[0])
        headingVectorY.append(heading[1])
        c, s = aspect[0], aspect[1]
        R = np.array(((c,s), (-s, c)))
        T = np.dot(R,np.transpose(heading))
        BetaX.append(T[0])
        BetaY.append(T[1])
        
    return Gradient, BetaX, BetaY

def ab2pitch(alpha,beta):
    return rad2deg*np.arccos(np.cos(deg2rad*alpha)/np.sqrt(np.cos(beta)**2+np.cos(deg2rad*alpha)**2*np.sin(beta)**2))


# =============================================================================
#    CAMIS cost function
# =============================================================================
def computeCAMIScost(B,Cd,Ca,Cl1,Cl2):
    Bt = np.transpose(B)
    K1 = (Ca+Cd)/2
    K2 = (Cl1+Cl2)/2
    D1 = (Ca-Cd)/2
    D2 = (Cl2-Cl1)/2
    Q1 = np.power(K1,2)
    Q2 = np.power(K2,2)
    Q3 = D1*D2
    Q = np.array(((Q1,Q3), (Q3, Q2)))
    D = np.array(((D1), (D2)))
    return (np.sqrt(np.dot(np.dot(Bt,Q),B))-np.dot(Bt,D))

def getCAMIScost(B,Q1,Q2,D1,D2):
    return np.sqrt(Q1*B[0]**2+Q2*B[1]**2+2*D1*D2*B[0]*B[1]) - \
           (B[0]*D1 + B[1]*D2)

# =============================================================================
#    Explicit formulation of CAMIS
# =============================================================================
           
def computeDirCosts(gradient, beta, cost, sigma):
    popt,_ = curve_fit(fittingCAMIS, (gradient,beta), cost, sigma=sigma)
#    cdRoots = (popt[0],popt[1],popt[-1])
#    caRoots = (popt[2],popt[3],popt[-1])
#    cl1Roots = (popt[4],popt[5],popt[-1])
#    cl2Roots = (popt[4],popt[5],popt[-1])
    cdRoots = (popt[0],popt[-1])
    caRoots = (popt[1],popt[-1])
    cl1Roots = (popt[2],popt[-1])
    cl2Roots = (popt[2],popt[-1])
    return cdRoots, caRoots, cl1Roots, cl2Roots 

def fittingCAMIS(x,x1,x2,x3,x4,x5,x6,Co):
    alpha, beta = x
#    Cd = dirCost(alpha, [x1, x2, Co])
#    Ca = dirCost(alpha, [x3, x4, Co])
#    Cl1 = dirCost(alpha, [x5, x6, Co])
#    Cl2 = dirCost(alpha, [x5, x6, Co])
    Cd = dirCost(alpha, [x1, Co])
    Ca = dirCost(alpha, [x2, Co])
    Cl1 = dirCost(alpha, [x3, Co])
    Cl2 = dirCost(alpha, [x3, Co])
    cBeta = np.cos(beta)
    sBeta = np.sin(beta)
    K1 = (Ca+Cd)/2
    K2 = (Cl1+Cl2)/2
    D1 = (Ca-Cd)/2
    D2 = (Cl2-Cl1)/2
    Q1 = np.power(K1,2)
    Q2 = np.power(K2,2)
    return np.sqrt(Q1*cBeta**2+Q2*sBeta**2+2*D1*D2*cBeta*sBeta) - \
           (cBeta*D1 + sBeta*D2)
    
def dirCost(gradient, K):
    # Gradient in degrees!
    return np.polyval(K, gradient)


def computeDrivingDirCosts(gradient, beta, cost, sigma):
#    bounds = ([0.0,0.0,0.0,0.0,0.0],[np.inf,1.0,np.inf,1.0,np.inf])
    bounds = ([0.01,0.01,0.01],[np.inf,np.inf,np.inf])
    p0 = ([20,20,20])
    popt,_ = curve_fit(fittingDrivingCAMIS, (gradient,beta), cost, sigma=sigma, bounds = bounds,p0 = p0,method='dogbox')
    # x1 = friction
    # x2 =
    return popt 

def computeAniCoLUT(cdRoots, caRoots, cl1Roots, cl2Roots, maxGradient):
    linearGradient = np.linspace(0,maxGradient,maxGradient+1)
    heading = np.arange(0, 2*np.pi, 0.01)
    aspect = [1,0] 
    Cd = dirCost(linearGradient, cdRoots)
    Ca = dirCost(linearGradient, caRoots)
    Cl1 = dirCost(linearGradient, cl1Roots)
    Cl2 = dirCost(linearGradient, cl2Roots)
    
    Cs = []
    anicoLUT = np.zeros((2,linearGradient.size))
    Anisotropy = np.zeros_like(linearGradient)
    for i,g in enumerate(linearGradient):
        for theta in heading:
            B = computeBeta(aspect,theta)
            preCost = computeCAMIScost(B,Cd[i],Ca[i],Cl1[i],Cl2[i])
            Cs.append(preCost)
        Anisotropy[i] = max(Cs)/min(Cs)
        Cs = []
    
    anicoLUT[:][0] = linearGradient
    anicoLUT[:][1] = Anisotropy
    
    return anicoLUT



# =============================================================================
#    Computation of the beta vector
# =============================================================================
def computeBeta(aspect,heading):
#    c, s = np.cos(aspect), np.sin(aspect)
#    R = np.array(((c,s), (-s, c)))
    if np.asarray(aspect).size == 1:
        R = np.array(((np.cos(aspect),np.sin(aspect)), (-np.sin(aspect), np.cos(aspect))))
    else:
        R = np.array(((aspect[0],aspect[1]), (-aspect[1], aspect[0])))
    if heading.size == 1:
        aHeading = np.array(((np.cos(heading)), (np.sin(heading))))
    else:
        aHeading = heading
    return np.dot(R,np.transpose(aHeading))

def getVectorialCostMap(slopeMap,CdRoots,CaRoots,Cl1Roots,Cl2Roots,AniCoLUT):
    vectorialCostMap = np.zeros([5,slopeMap.shape[0],slopeMap.shape[1]])
    for i in range(slopeMap.shape[1]):
        for j in range(slopeMap.shape[0]):
            if (slopeMap[j][i] > AniCoLUT[0][-1]):
#                Cobs = dirCost(AniCoLUT[0][-1], CdRoots) + 2*(slopeMap[j][i]-AniCoLUT[0][-1])**2
                Cobs = np.max((dirCost(slopeMap[j][i], CaRoots),dirCost(slopeMap[j][i], Cl1Roots),dirCost(slopeMap[j][i], Cl2Roots)))
                vectorialCostMap[0][j][i] = 1.0
                vectorialCostMap[1][j][i] = Cobs**2 # Q1
                vectorialCostMap[2][j][i] = Cobs**2# Q2
                vectorialCostMap[3][j][i] = 0.0 # D1
                vectorialCostMap[4][j][i] = 0.0 # D2
            else:
                Cd = dirCost(slopeMap[j][i], CdRoots)
                Ca = dirCost(slopeMap[j][i], CaRoots)
                Cl1 = dirCost(slopeMap[j][i], Cl1Roots)
                Cl2 = dirCost(slopeMap[j][i], Cl2Roots)
                vectorialCostMap[0][j][i] = getAnisotropy(slopeMap[j][i],AniCoLUT)
                vectorialCostMap[1][j][i] = ((Ca+Cd)/2)**2 # Q1
                vectorialCostMap[2][j][i] = ((Cl1+Cl2)/2)**2# Q2
                vectorialCostMap[3][j][i] =  (Ca-Cd)/2 # D1
                vectorialCostMap[4][j][i] = (Cl2-Cl1)/2 # D2    
    return vectorialCostMap

def getAnisotropy(slope,AniCoLUT):
    if slope <= AniCoLUT[0][0]:
        return AniCoLUT[1][0]
    if slope >= AniCoLUT[0][-1]:
        return AniCoLUT[1][-1]
    index = bisect_left(AniCoLUT[0][:], slope)
    decimalPart = (slope - AniCoLUT[0][index-1])/(AniCoLUT[0][index] - AniCoLUT[0][index-1])
    anisotropy = decimalPart*(AniCoLUT[1][index]-AniCoLUT[1][index-1]) + AniCoLUT[1][index-1]
    return anisotropy


def saveCamis(fileName,CdRoots,CaRoots,Cl1Roots,Cl2Roots,AniCoLUT):
    with open(fileName, mode='w') as file:
        camis_writer = csv.writer(file, delimiter=',', quotechar='"',\
                                  quoting=csv.QUOTE_MINIMAL)
        camis_writer.writerow(CdRoots)
        camis_writer.writerow(CaRoots)
        camis_writer.writerow(Cl1Roots)
        camis_writer.writerow(Cl2Roots)
        camis_writer.writerow(AniCoLUT[0][:])
        camis_writer.writerow(AniCoLUT[1][:])
    return

def readCamis(fileName):
    with open(fileName) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',', quotechar='"',\
                                  quoting=csv.QUOTE_MINIMAL)
        CdRoots = [float(i) for i in next(csv_reader)]
        next(csv_reader)
        CaRoots = [float(i) for i in next(csv_reader)]
        next(csv_reader)
        Cl1Roots = [float(i) for i in next(csv_reader)]
        next(csv_reader)
        Cl2Roots = [float(i) for i in next(csv_reader)]
        next(csv_reader)
        gradientArray = [float(i) for i in next(csv_reader)]
        next(csv_reader)
        aniCoeff = [float(i) for i in next(csv_reader)]
        AniCoLUT = np.zeros((2,len(gradientArray)))
        AniCoLUT[:][0] = gradientArray
        AniCoLUT[:][1] = aniCoeff
    return CdRoots,CaRoots,Cl1Roots,Cl2Roots,AniCoLUT

def showCAMIS(CdRoots, CaRoots, Cl1Roots, Cl2Roots,beta,gradient,cost):
    
#    Cd = asymptoticCost(linearGradient,27.0, CdRoots[0], CdRoots[-1])
#    Ca = asymptoticCost(linearGradient, 27.0, CaRoots[0], CaRoots[-1])
#    Cl1 = asymptoticCost(linearGradient, 27.0, Cl1Roots[0], Cl1Roots[-1])
#    Cl2 = asymptoticCost(linearGradient, 27.0, Cl2Roots[0], Cl2Roots[-1])
    
    linearGradient = np.linspace(0,30,31)
    heading = np.arange(0, 2*np.pi, 0.01)
    aspect = [1,0]
    
    
    Cd = dirCost(linearGradient, CdRoots)
    Ca = dirCost(linearGradient, CaRoots)
    Cl1 = dirCost(linearGradient, Cl1Roots)
    Cl2 = dirCost(linearGradient, Cl2Roots)
    
    CMax = max(Ca[-2],Cl1[-2])
    CMax = max(CMax,Cl2[-2])
    costX = np.multiply(np.cos(beta),cost)
    costY = np.multiply(np.sin(beta),cost)
    
    rgba_colors = np.zeros((costX.shape[0],4))
    rgba_colors[:,0] = np.abs(np.sin(beta))
    rgba_colors[:, 3] = np.abs(np.sin(beta))
    
    
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
    
    fig1 = plt.figure(figsize=(8, 8))
    ax1 = plt.axes(rect_scatter)
    ax1.tick_params(direction='in', top=True, right=True)
    ax1.set_facecolor('xkcd:pale grey')
    ax1.set_title('$\longleftarrow$ Lateral Direction #2      Lateral Direction #1 $\longrightarrow$',y = -0.13)
    plt.grid('True')
    ax1.fill_betweenx(linearGradient,0,Cl1,color='xkcd:sea blue')
    ax1.fill_betweenx(linearGradient,0,-Cl2,color='xkcd:leaf')
    ax1.scatter(costY,gradient,color=rgba_colors)
    ax1.legend(['Lateral Cost 1','Lateral Cost 2','Experimental Data'])
    
    ax1.set_xlim(-CMax,CMax)
    gradientArray = np.asarray(gradient)
    ax1.set_ylim(0,gradientArray.max())
    ax1.plot(Cl1,linearGradient, lw = 2, color = 'xkcd:dusk blue')
    ax1.plot(-Cl2,linearGradient, lw = 2, color = 'xkcd:camo green')
    plt.xlabel('Cost [As/m]')
    plt.ylabel('Slope Gradient [degrees]')
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)
    
    binwidth = 0.25
    binsX = np.arange(-CMax, CMax + binwidth, binwidth)
    binsY = np.arange(0, gradientArray.max() + binwidth, binwidth)
    ax_histx.hist(costY, bins=binsX, color = "m")
    ax_histy.hist(gradient, bins=binsY, orientation='horizontal', color = "m")
    
    ax_histx.set_xlim(ax1.get_xlim())
    ax_histy.set_ylim(ax1.get_ylim())


    fig1.tight_layout()
    
    rgba_colors[:,0] = np.abs(np.cos(beta))
    rgba_colors[:, 3] = np.abs(np.cos(beta))
    
    
    fig2 = plt.figure(figsize=(8, 8))
    ax2 = plt.axes(rect_scatter)
    ax2.set_facecolor('xkcd:pale grey')
    ax2.set_title('$\longleftarrow$ Ascent Direction      Descent Direction $\longrightarrow$',y = -0.13)
    plt.grid('True')
    ax2.fill_betweenx(linearGradient,0,Cd,color='xkcd:sea blue')
    ax2.fill_betweenx(linearGradient,0,-Ca,color='xkcd:leaf')
    ax2.scatter(costX,gradient,color=rgba_colors)
    ax2.legend(['Descent Cost','Ascent Cost','Experimental Data'])
    ax2.set_xlim(-CMax,CMax)
    ax2.set_ylim(0,gradientArray.max())
    ax2.plot(Cd,linearGradient, lw = 2, color = 'xkcd:dusk blue')
    ax2.plot(-Ca,linearGradient, lw = 2, color = 'xkcd:camo green')
    plt.xlabel('Cost [As/m]')
    plt.ylabel('Slope Gradient [degrees]')
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)
    
    binwidth = 0.25
    binsX = np.arange(-CMax, CMax + binwidth, binwidth)
    binsY = np.arange(0, gradientArray.max() + binwidth, binwidth)
    ax_histx.hist(costX, bins=binsX, color = "m")
    ax_histy.hist(gradient, bins=binsY, orientation='horizontal', color = "m")
    
    ax_histx.set_xlim(ax2.get_xlim())
    ax_histy.set_ylim(ax2.get_ylim())
    
    fig2.tight_layout()
    
    
    
    
    
#    fig3 = plt.figure()
#    
#    sqrtQ1 = np.zeros_like(linearGradient)
#    sqrtQ2 = np.zeros_like(linearGradient)
#    D1 = np.zeros_like(linearGradient)
#    D2 = np.zeros_like(linearGradient)
#    for i,g in enumerate(linearGradient):
#        sqrtQ1[i] = (Ca[i]+Cd[i])/2
#        sqrtQ2[i] = (Cl2[i]+Cl1[i])/2
#        D1[i] = (Ca[i]-Cd[i])/2
#        D2[i] = (Cl2[i]-Cl1[i])/2
#    ax3 = fig3.add_subplot(1, 1, 1)
#    ax3.plot(linearGradient,sqrtQ1,color='b')
#    ax3.plot(linearGradient,sqrtQ2,color='g')
#    ax3.plot(linearGradient,D1,color='m')
#    ax3.plot(linearGradient,D2,color='y')
#    ax3.plot(linearGradient,Ca,color='burlywood')
#    ax3.plot(linearGradient,Cd,color='lightblue')
#    ax3.plot(linearGradient,Cl1,color='yellowgreen')
#    ax3.plot(linearGradient,Cl2,color='coral')
#    ax3.legend(['sqrt(Q1)','sqrt(Q2)','D1','D2','Ascent Cost Ca',\
#                'Descent Cost Cd','Lateral Cost Cl1','Lateral Cost Cl2'])
#    plt.ylabel('Cost [As/m]')
#    plt.xlabel('Slope Gradient [degrees]')
#    
#    Bs = []
#    Cs = []
#    Anisotropy = np.zeros_like(linearGradient)
#    for i,g in enumerate(linearGradient):
#        for theta in heading:
#            B = computeBeta(aspect,theta)
#            Bs.append(np.arctan2(B[1],B[0]))
#            preCost = computeCAMIScost(B,Cd[i],Ca[i],Cl1[i],Cl2[i])
#            Cs.append(preCost)
#        Anisotropy[i] = max(Cs)/min(Cs)
#        Cs = []
#        Bs = []
#    ax3b = ax3.twinx()
#    ax3b.plot(linearGradient,Anisotropy,color='r')
#    ax3b.set_ylabel('Anisotropy', color='r')
#    ax3b.tick_params('y', colors='r')
#    fig3.tight_layout()
#    
#
#
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    ax.plot(Cd, np.zeros_like(linearGradient), linearGradient, color = 'r')
#    ax.plot(-Ca, np.zeros_like(linearGradient), linearGradient, color = 'r')
#    ax.plot(np.zeros_like(linearGradient), Cl1, linearGradient, color = 'r')
#    ax.plot(np.zeros_like(linearGradient), -Cl2, linearGradient, color = 'r')
#    Xs = []
#    Ys = []
#     
#    for i,g in enumerate(linearGradient):
#        for theta in heading:
#            B = computeBeta(aspect,theta)
#            preCost = computeCAMIScost(B,Cd[i],Ca[i],Cl1[i],Cl2[i])
#            Xs.append(B[0]*preCost)
#            Ys.append(B[1]*preCost)
#        ax.plot(Xs, Ys, g, color=plt.cm.jet(float(g)/25))
#        Xs = []
#        Ys = []
#    ax.set_xlim(-CMax,CMax)
#    ax.set_ylim(-CMax,CMax)
#    ax.set_aspect('equal')
#    fig.tight_layout()
#    
#    fig5 = plt.figure()
#    ax5 = fig5.add_subplot(1, 1, 1)
#    ax5 = plt.subplot(111, projection='polar')
#    ax5.set_facecolor('xkcd:light blue')
#    Bs = []
#    Cs = []
#    for i,g in enumerate(linearGradient):
#        for theta in heading:
#            B = computeBeta(aspect,theta)
#            Bs.append(np.arctan2(B[1],B[0]))
#            preCost = computeCAMIScost(B,Cd[i],Ca[i],Cl1[i],Cl2[i])
#            Cs.append(preCost)
#        ax5.plot(Bs, Cs, 'xkcd:sea blue', lw = 2)
#        Cs = []
#        Bs = []
#    fig5.tight_layout()
        
    
    
    
    
    
    
    
    
        
#    xs = np.multiply(np.cos(beta),cost)
#    ys = np.multiply(np.sin(beta),cost)
#    col = np.zeros_like(xs)
#    for i,s in enumerate(xs):
#        col[i] = np.sqrt(np.power(xs[i],2)+np.power(ys[i],2))
#    ax.scatter(xs, ys, gradient,c = col, alpha = .1)
    
    
    plt.show()
        
    
    
# =============================================================================
#    Class of CAMIS model for any robot
# =============================================================================
    
class CamisDrivingModel:
    def __init__(self, robot_data):
        self.slopeThreshold = robot_data['slope_threshold']
        self.occupancy_radius = robot_data['occupancy_radius']
        self.tracking_error = robot_data['tracking_error']
        self.speed = robot_data['speed']
        self.mode = robot_data['mode']
        self.rolling_resistance_mode = robot_data['rolling_resistance_mode']
        self.friction = robot_data['friction']
        self.kmg = robot_data['kmg']
        self.steepness_brake_margin = robot_data['steepness_brake_margin']
        self.risk_mode = robot_data['risk_mode']
        self.risk_roll = robot_data['risk_roll']*deg2rad
        self.risk_pitch = robot_data['risk_pitch']*deg2rad
        self.risk_block = robot_data['risk_block']*deg2rad
        self.risk_init = robot_data['risk_init']*deg2rad
        self.risk_mult = robot_data['risk_mult']
        self.risk_margin = robot_data['risk_margin']*deg2rad
        self.bezier_coeff = np.minimum(np.maximum(robot_data['bezier_coeff'], 0.0), 0.99)
        self.slip_mode = robot_data['slip_mode']
        if (self.slip_mode == 'roots'):
            self.slip_roots_ascent = robot_data['slip_roots_ascent']
            self.slip_roots_descent = robot_data['slip_roots_descent']
            self.slip_roots_lateral = robot_data['slip_roots_lateral']
        self.slip = robot_data['slip']
        self.roll_weight = robot_data['roll_weight']
        self.ascent_weight = robot_data['ascent_weight']
        self.descent_weight = robot_data['descent_weight']
        self.computeBezierPoints()
        self.computeAniCoLUT()
    
    def setAsIsotropic(self):
        self.mode = 'isotropic'
        self.computeBezierPoints()
        self.computeAniCoLUT()
        
    def fitCAMIS(self, gradient,beta,cost,sigma):
        bounds = ([0.01,0.01,0.01,0.0,0.0],[1.0,1.0,np.inf,self.slopeThreshold/45.0,self.slopeThreshold/45.0])
        popt,_ = curve_fit(self.fittingDrivingCAMIS, (gradient,beta), cost, sigma=sigma, bounds = bounds,method='dogbox')
        self.friction_parallel = popt[0]
        self.friction_perp = popt[1]
        self.kmg = popt[2]
        self.slip_parallel = popt[3]
        self.slip_perp = popt[4]
#        self.fmg_parallel = popt[0]
#        self.fmg_perp = popt[1]
#        self.mg = popt[2]
#        self.friction_parallel = self.fmg_parallel/self.mg
#        self.friction_perp = self.fmg_perp/self.mg
        self.computeAniCoLUT() 
        
    def fittingDrivingCAMIS(self,x,x1,x2,x3,x4,x5):
        alpha, beta = x
        cBeta = np.cos(beta)
        sBeta = np.sin(beta)
        K1 = x1*x3*np.cos(deg2rad*alpha)/(1-np.sin(deg2rad*alpha*x4))
        K2 = x2*x3*np.cos(deg2rad*alpha)/(1-np.sin(deg2rad*alpha*x4))
        D1 = x3*np.sin(deg2rad*alpha)/(1-np.sin(deg2rad*alpha*x5))
        D2 = x3*np.sin(deg2rad*alpha)/(1-np.sin(deg2rad*alpha*x5))
        Q1 = np.power(K1,2)
        Q2 = np.power(K2,2)
        return np.sqrt(Q1*cBeta**2+Q2*sBeta**2+2*D1*D2*cBeta*sBeta) - \
               (cBeta*D1 + sBeta*D2)
           
    def getVectorialCostMap(self,slopeMap):
        vectorialCostMap = np.ones([5,slopeMap.shape[0],slopeMap.shape[1]])
        Cobs = np.max((self.getCd(self.slopeThreshold),\
                                   self.getCa(self.slopeThreshold),\
                                   self.getCl(self.slopeThreshold)))*self.getAnisotropy(self.slopeThreshold)
        vectorialCostMap[1] = vectorialCostMap[1]*Cobs**2
        vectorialCostMap[2] = vectorialCostMap[2]*Cobs**2
        vectorialCostMap[3] = vectorialCostMap[3]*0.0
        vectorialCostMap[4] = vectorialCostMap[4]*0.0
        
        for i in range(slopeMap.shape[1]):
            for j in range(slopeMap.shape[0]):
                slope = slopeMap[j][i]
#                if (slope < self.slopeThreshold):
                if (self.mode == 'isotropic'):
                        Cn = self.getCn(slope)
                        vectorialCostMap[0][j][i] = 1
                        vectorialCostMap[1][j][i] = Cn**2 # Q1
                        vectorialCostMap[2][j][i] = Cn**2# Q2
                        vectorialCostMap[3][j][i] =  0
                        vectorialCostMap[4][j][i] = 0
                else:
                        Cd = self.getCd(slope)
                        Ca = self.getCa(slope)
                        Cl1 = self.getCl(slope)
                        Cl2 = self.getCl(slope)
                        vectorialCostMap[0][j][i] = self.getAnisotropy(slope)
                        vectorialCostMap[1][j][i] = ((Ca+Cd)/2)**2 # Q1
                        vectorialCostMap[2][j][i] = ((Cl1+Cl2)/2)**2# Q2
                        vectorialCostMap[3][j][i] =  (Ca-Cd)/2 # D1
                        vectorialCostMap[4][j][i] = (Cl2-Cl1)/2 # D2    
        return vectorialCostMap
    
    def getQuadraticBezierCost(self,x, P0, P1, P2):
        X1 = P0[0] - P1[0]
        X2 = P2[0] - P1[0]
        Y1 = P0[1] - P1[1]
        Y2 = P2[1] - P1[1]
        
        if x < P0[0]:
            return P0[1] - (P0[0] - x)/(P1[0] - P0[0])*(P1[1] - P0[1])
        
        if x > P2[0]:
            return P2[1] + (x - P2[0])/(P2[0] - P1[0])*(P2[1] - P1[1])
        else:
            if np.abs(X1 - X2) < 0.001:
                t = (x - P0[0] - P1[0] + X1)/(2*X1)
            elif np.abs(-X1 - X2) < 0.001:
                t = (X1 + P1[0] - x)/(2*X1)
            else:
                t = (-2*X1 - np.sqrt(4*X1**2 - 
                                     4*(x-P1[0]-X1)*(-X1-X2)))/(2*(-X1-X2))
            Cost = P1[1] + (1 - t)**2*Y1 + t**2*Y2
            return Cost
    def getDescentBezier(self,steepness):
        t = (self.brakePoint01[0] - steepness)/(-2*self.steepness_brake_margin*deg2rad)
        return (1 - t)**2*self.brakePoint01[1] + t**2*self.brakePoint02[1]
    
    def getCubicBezierCost(self,steepness,P0,P1,P2,P3):
#        res = scipy.optimize.least_squares(lambda x : 
#        (1-x)**3*P0[0] + 3*(1-x)**2*x*P1[0] + 3*(1-x)*x**2*P2[0]+x**3*P3[0] - 
#            steepness,0.5, bounds = (0,1))
        coeff = [-P0[0]+3*P1[0]-3*P2[0]+P3[0], 
                 3*P0[0] - 6*P1[0] + 3*P2[0],
                 -3*P0[0] + 3*P1[0],
                 P0[0] - steepness]
        results = np.roots(coeff)
        for res in results:
            if (res >= 0.0)and(res <=1.0)and(np.isreal(res)):
                t = res.real
#        t = res.x
        return np.asscalar((1-t)**3*P0[1] + 3*(1-t)**2*t*P1[1] + 3*(1-t)*t**2*P2[1] + t**3*P3[1])
                
    def computeBezierPoints(self):
        # The point at which the robot starts braking
        self.brakeSteepness_rad = np.arctan(self.friction)
        self.brakePoint = np.array([self.brakeSteepness_rad,0])
        
        # Computing the Bezier Braking points
        if self.brakeSteepness_rad < self.steepness_brake_margin*deg2rad:
            brake01X = 0
        else:
            brake01X = (self.brakeSteepness_rad - self.steepness_brake_margin*deg2rad)
        brake01Y = self.kmg * (self.friction - np.tan(brake01X))
        self.brakePoint01 = np.array([brake01X, brake01Y])
        
        brake02X = (self.brakeSteepness_rad + self.steepness_brake_margin*deg2rad)
        brake02Y = self.kmg * np.abs(self.friction - np.tan(brake02X))
        self.brakePoint02 = np.array([brake02X, brake02Y]) 
        
        
        #Get the risk line
        self.ascentRiskThreshold = np.min([self.risk_pitch,self.risk_block])
        self.Cmax = self.risk_mult*self.getRawCa(self.ascentRiskThreshold*rad2deg)
        self.risk_gain = self.Cmax/(self.bezier_coeff*self.risk_margin)
        
        #Obtain the ascent Bezier points
        initialPointX = self.risk_init
        initialPointY = self.getRawCa(initialPointX*rad2deg)
        
        aa = self.ascentRiskThreshold - self.risk_init
        
        intersection1X = self.risk_init + self.bezier_coeff*aa/2
        intersection1Y = self.getRawCa(intersection1X*rad2deg)

        intersection2X = self.ascentRiskThreshold - self.bezier_coeff*aa/2
        intersection2Y = self.Cmax
        
        self.ascentBezierPoint_initial = np.array([initialPointX,initialPointY])
        self.ascentBezierPoint_intersection1 = np.array([intersection1X,intersection1Y])
        self.ascentBezierPoint_intersection2 = np.array([intersection2X,intersection2Y])
        self.ascentBezierPoint_risk = np.array([self.ascentRiskThreshold,self.Cmax])
        
        #Obtain the lateral Bezier points
        self.lateralRiskThreshold = self.risk_roll
        initialPointX = self.risk_init
        initialPointY = self.getRawCl(initialPointX*rad2deg)
        
        aa = self.lateralRiskThreshold - self.risk_init
        intersection1X = self.risk_init + self.bezier_coeff*aa/2
        intersection1Y = self.getRawCl(intersection1X*rad2deg)
        
        intersection2X = self.lateralRiskThreshold - self.bezier_coeff*aa/2
        intersection2Y = self.Cmax
        
        self.lateralBezierPoint_initial = np.array([initialPointX,initialPointY])
        self.lateralBezierPoint_intersection1 = np.array([intersection1X,intersection1Y])
        self.lateralBezierPoint_intersection2 = np.array([intersection2X,intersection2Y])
        self.lateralBezierPoint_risk = np.array([self.lateralRiskThreshold,self.Cmax])
        
        #Obtain the descent Bezier points
        self.descentRiskThreshold = self.risk_pitch
        initialPointX = self.risk_init
        initialPointY = self.getRawCd(initialPointX*rad2deg)
        
        aa = self.descentRiskThreshold - self.risk_init
        intersection1X = self.risk_init + self.bezier_coeff*aa/2
        intersection1Y = self.getRawCd(intersection1X*rad2deg)
        
        intersection2X = self.descentRiskThreshold - self.bezier_coeff*aa/2
        intersection2Y = self.Cmax
        
        self.descentBezierPoint_initial = np.array([initialPointX,initialPointY])
        self.descentBezierPoint_intersection1 = np.array([intersection1X,intersection1Y])
        self.descentBezierPoint_intersection2 = np.array([intersection2X,intersection2Y])
        self.descentBezierPoint_risk = np.array([self.descentRiskThreshold,self.Cmax])
        
        #The intersection between risk function and ascent function
        intersectionX = scipy.optimize.fsolve(lambda x : 
        (x - np.minimum(self.risk_pitch, self.risk_block) + 
         self.risk_margin)/self.risk_margin*self.risk_gain - np.tan(x) - 
         self.friction,0)
        intersectionY = self.friction + np.tan(intersectionX)
        self.ascentIntersection = np.array([intersectionX,intersectionY])
        
        #The intersection between risk function and descent function
        if self.brakeSteepness_rad >= self.risk_pitch - self.risk_margin:
            intersection2X = scipy.optimize.fsolve(lambda x : 
                (x - self.risk_pitch + 
                 self.risk_margin)/self.risk_margin*self.risk_gain + 
                 np.tan(x) - self.friction, 0)
            intersection2Y = self.friction - np.tan(intersection2X)
        else:
            intersection2X = scipy.optimize.fsolve(lambda x : 
                (x - self.risk_pitch + 
                 self.risk_margin)/self.risk_margin*self.risk_gain + 
                -np.tan(x) + self.friction,0)
            intersection2Y = - self.friction + np.tan(intersection2X)
        self.descentIntersection = np.array([intersection2X,intersection2Y])
        
        #The intersection between risk function and lateral function 
        intersection3X = scipy.optimize.fsolve(lambda x : 
                (x - self.risk_roll + 
                 self.risk_margin)/self.risk_margin*self.risk_gain - 
                 self.friction*np.cos(x), 0)
        intersection3Y = np.sqrt(self.friction**2  - 
                                 np.tan(intersection3X)**2)

        self.lateralIntersection = np.array([intersection3X,intersection3Y])
        
        #The first Bezier points
        self.initialPoint = np.array([(1-self.bezier_coeff)*intersectionX,
                                 self.friction + 
                                 np.tan((1-self.bezier_coeff)*intersectionX)])
        
        if self.brakeSteepness_rad >= self.risk_pitch - self.risk_margin:
            self.initial2Point = np.array([(1-self.bezier_coeff)*intersection2X,
                                      self.friction - 
                                      np.tan((1-self.bezier_coeff)*intersection2X)])
        else:
            self.initial2Point = np.array([(1-self.bezier_coeff)*self.brakeSteepness_rad,
                                      self.friction - np.tan((1-self.bezier_coeff)*self.brakeSteepness_rad)])
            
        self.initial3Point = np.array([(1-self.bezier_coeff)*intersection3X,
                                       self.friction*np.cos((1-self.bezier_coeff)*intersection3X)])
            
        # The last Bezier points
        self.riskPoint = np.array([intersectionX + self.bezier_coeff*(np.minimum(self.risk_block, self.risk_pitch) - intersectionX),
                              intersectionY + self.bezier_coeff*(self.risk_gain-intersectionY)])
        self.riskPoint2 = np.array([intersection2X + self.bezier_coeff*(self.risk_pitch - intersection2X),
                               intersection2Y + self.bezier_coeff*(self.risk_gain-intersection2Y)])
        self.riskPoint3 = np.array([intersection3X + self.bezier_coeff*(self.risk_roll - intersection3X),
                               intersection3Y + self.bezier_coeff*(self.risk_gain-intersection3Y)])

#    def getSlipFactor(self, steepness):
#        if self.slip_mode == 'none':
#            return 1
#        elif self.slip_mode == 'sin':
#            return 1/(1 - np.sin(np.min([steepness/self.slopeThreshold,0.99])*90.0*deg2rad)**self.slip)
#        elif self.slip_mode == 'roots':
#            s = np.min([np.polyval(self.slip_roots,steepness),0.95])
#            return 1/(1 - s)
#        else:
#            raise ValueError('The slip mode ' + self.slip_mode + ' is not valid')
    
    def getSSa(self, steepness_deg):
        if self.slip_mode == 'none':
            return 1
        elif self.slip_mode == 'sin':
            return 1/(1 - np.sin(np.min([steepness_deg/self.slopeThreshold,0.99])*90.0*deg2rad)**self.slip)
        elif self.slip_mode == 'roots':
            if steepness_deg >= self.slopeThreshold:
                s = np.max([np.min([np.polyval(self.slip_roots_ascent,self.slopeThreshold),0.99]),
                            np.min([np.polyval(self.slip_roots_lateral,self.slopeThreshold),0.99]),
                            np.min([np.polyval(self.slip_roots_descent,self.slopeThreshold),0.99])])
            else:
                s = np.min([np.polyval(self.slip_roots_ascent,steepness_deg),0.99])
            return 1/(1 - s)
        else:
            raise ValueError('The slip mode ' + self.slip_mode + ' is not valid')
            
    def getSSl(self, steepness_deg):
        if self.slip_mode == 'none':
            return 1
        elif self.slip_mode == 'sin':
            return 1/(1 - np.sin(np.min([steepness_deg/self.slopeThreshold,0.99])*90.0*deg2rad)**self.slip)
        elif self.slip_mode == 'roots':
            if steepness_deg >= self.slopeThreshold:
                s = np.max([np.min([np.polyval(self.slip_roots_ascent,self.slopeThreshold),0.99]),
                            np.min([np.polyval(self.slip_roots_lateral,self.slopeThreshold),0.99]),
                            np.min([np.polyval(self.slip_roots_descent,self.slopeThreshold),0.99])])
            else:
                s = np.min([np.polyval(self.slip_roots_lateral,steepness_deg),0.99])
            return 1/(1 - s)
        else:
            raise ValueError('The slip mode ' + self.slip_mode + ' is not valid')
            
    def getSSd(self, steepness_deg):
        if self.slip_mode == 'none':
            return 1
        elif self.slip_mode == 'sin':
            return 1/(1 - np.sin(np.min([steepness_deg/self.slopeThreshold,0.99])*90.0*deg2rad)**self.slip)
        elif self.slip_mode == 'roots':
            if steepness_deg >= self.slopeThreshold:
                s = np.max([np.min([np.polyval(self.slip_roots_ascent,self.slopeThreshold),0.99]),
                            np.min([np.polyval(self.slip_roots_lateral,self.slopeThreshold),0.99]),
                            np.min([np.polyval(self.slip_roots_descent,self.slopeThreshold),0.99])])
            else:
                s = np.min([np.polyval(self.slip_roots_descent,steepness_deg),0.99])
            return 1/(1 - s)
        else:
            raise ValueError('The slip mode ' + self.slip_mode + ' is not valid')
            
    def getRRd(self, steepness_deg):
        steepness = steepness_deg*deg2rad
        if steepness < 0:
            raise ValueError('ERROR: input value of steepness is negative')
        if steepness < self.brakePoint01[0]:
            return (self.friction - np.tan(steepness))*self.kmg
        elif steepness > self.brakePoint02[0]:
            return (- self.friction + np.tan(steepness))*self.kmg
        else:
            return self.getDescentBezier(steepness)
#            return self.getQuadraticBezierCost(steepness, self.brakePoint01, 
#                                               self.brakePoint, self.brakePoint02)
        
    def getRRa(self, steepness_deg):
        steepness = steepness_deg*deg2rad
        if steepness < 0:
            raise ValueError('ERROR: input value of steepness is negative')
        else:
            return (self.friction + np.tan(steepness))*self.kmg
        
    def getRRl(self, steepness_deg):
        steepness = steepness_deg*deg2rad
        if steepness < 0:
            raise ValueError('ERROR: input value of steepness is negative')
        else:
            return (self.friction * np.cos(steepness))*self.kmg
#            return np.sqrt(self.getRRa(steepness_deg)*self.getRRd(steepness_deg))
#            return (self.getRRa(steepness_deg) + self.getRRd(steepness_deg))/2
        
    def getRawCd(self, steepness_deg):
        S = self.getSSd(steepness_deg)
        if self.rolling_resistance_mode == 'none':
            R = 1
        else:
            R = self.getRRd(steepness_deg)
        pv = self.speed*np.cos(steepness_deg*deg2rad)
        return R * S / pv
#        return 0.1*self.getRawCa(steepness_deg) + 0.9*self.getRawCl(steepness_deg)
    def getCd(self, steepness_deg):
        rawC = self.getRawCd(steepness_deg)
        W = (1 + self.descent_weight*np.tan(steepness_deg*deg2rad))
        return rawC * W
#        if self.risk_mode == 'none':
#            return rawC
#        if steepness_deg < self.descentBezierPoint_initial[0]*rad2deg:
#            return rawC
#        elif steepness_deg < self.descentBezierPoint_risk[0]*rad2deg:
#            return self.getQuadraticBezierCost(steepness_deg*deg2rad, 
#                                               self.descentBezierPoint_initial,
#                                               self.descentBezierPoint_intersection1,
#                                               self.descentBezierPoint_risk)
##            return self.getCubicBezierCost(steepness_deg*deg2rad, 
##                                               self.descentBezierPoint_initial,
##                                               self.descentBezierPoint_intersection1,
##                                               self.descentBezierPoint_intersection2,
##                                               self.descentBezierPoint_risk)
#        else:
#            return self.risk_gain*(steepness_deg*deg2rad - self.descentRiskThreshold) + self.Cmax
#            return self.Cmax
    
    def getRawCa(self,steepness_deg):
        S = self.getSSa(steepness_deg)
        if self.rolling_resistance_mode == 'none':
            R = 1
        else:
            R = self.getRRa(steepness_deg)
        pv = self.speed*np.cos(steepness_deg*deg2rad)
        return R * S / pv
    def getCa(self, steepness_deg):
        rawC = self.getRawCa(steepness_deg)
        W = (1 + self.ascent_weight*np.tan(steepness_deg*deg2rad))
        return rawC * W
#        if self.risk_mode == 'none':
#            return rawC
#        if steepness_deg < self.ascentBezierPoint_initial[0]*rad2deg:
#            return rawC
#        elif steepness_deg < self.ascentBezierPoint_risk[0]*rad2deg:
#            return self.getQuadraticBezierCost(steepness_deg*deg2rad, 
#                                               self.ascentBezierPoint_initial,
#                                               self.ascentBezierPoint_intersection1,
#                                               self.ascentBezierPoint_risk)
#        else:
#            return self.risk_gain*(steepness_deg*deg2rad - self.ascentRiskThreshold) + self.Cmax
#            return self.Cmax
        
    def getRawCl(self,steepness_deg):
        S = self.getSSl(steepness_deg)
        if self.rolling_resistance_mode == 'none':
            R = 1
        else:
            R = self.getRRl(steepness_deg)
        pv = self.speed
        return R * S / pv 
#        return 0.9*self.getRawCa(steepness_deg) + 0.1*R * S / pv
    def getCl(self, steepness_deg):
        rawC = self.getRawCl(steepness_deg)
        W = (1 + self.roll_weight*np.tan(steepness_deg*deg2rad))
        return rawC * W
#        if self.perp_coeff <= 1.0:
#            self.perp_coeff = np.max([self.perp_coeff,0.0])
#            return (1 - self.perp_coeff)*rawC + self.perp_coeff*self.getRawCa(steepness_deg)
#        else:
#            return self.getRawCa(steepness_deg)*(1 + (self.perp_coeff - 1)*np.tan(steepness_deg*deg2rad))
#        if self.risk_mode == 'none':
#            return self.getRawCa(steepness_deg)*(1 + np.tan(steepness_deg*deg2rad))
#        if steepness_deg < self.lateralBezierPoint_initial[0]*rad2deg:
#            return rawC
#        elif steepness_deg < self.lateralBezierPoint_risk[0]*rad2deg:
#            return self.getQuadraticBezierCost(steepness_deg*deg2rad, 
#                                               self.lateralBezierPoint_initial,
#                                               self.lateralBezierPoint_intersection1,
#                                               self.lateralBezierPoint_risk)
#        else:
#            return self.risk_gain*(steepness_deg*deg2rad - self.lateralRiskThreshold) + self.Cmax
#            return self.Cmax
        
    
    def getCn(self,slope): #TODO - check this
        cd = self.getCd(slope)
        ca = self.getCa(slope)
        cl = self.getCl(slope)
        return (ca + cd)/2*np.sqrt(cl/np.sqrt(ca*cd))
    
    def computeAniCoLUT(self):
        linearGradient = np.linspace(0,89,90)
        heading = np.arange(0, 2*np.pi, 0.01)
        aspect = [1,0] 
        
        Cs = []
        Cd = np.zeros_like(linearGradient)
        Ca = np.zeros_like(linearGradient)
        Cl1 = np.zeros_like(linearGradient)
        Cl2 = np.zeros_like(linearGradient)
        anicoLUT = np.zeros((2,linearGradient.size))
        Anisotropy = np.zeros_like(linearGradient)
        for i,g in enumerate(linearGradient):
            Cd[i] = self.getCd(linearGradient[i])
            Ca[i] = self.getCa(linearGradient[i])
            Cl1[i] = self.getCl(linearGradient[i])
            Cl2[i] = self.getCl(linearGradient[i])
            for theta in heading:
                B = computeBeta(aspect,theta)
                preCost = computeCAMIScost(B,Cd[i],Ca[i],Cl1[i],Cl2[i])
                Cs.append(preCost)
            Anisotropy[i] = max(Cs)/min(Cs)
            Cs = []
        
        anicoLUT[:][0] = linearGradient
        anicoLUT[:][1] = Anisotropy
    
        self.anicoLUT = anicoLUT
        
    def getAnisotropy(self,slope):
        if slope <= self.anicoLUT[0][0]:
            return self.anicoLUT[1][0]
        if slope >= self.anicoLUT[0][-1]:
            return self.anicoLUT[1][-1]
        index = bisect_left(self.anicoLUT[0][:], slope)
        decimalPart = (slope - self.anicoLUT[0][index-1])/(self.anicoLUT[0][index] - self.anicoLUT[0][index-1])
        anisotropy = decimalPart*(self.anicoLUT[1][index]-self.anicoLUT[1][index-1]) + self.anicoLUT[1][index-1]
        return anisotropy
    
    def getCost(self, slope, aspect, heading):
        beta = computeBeta(aspect,heading)
        Cd = self.getCd(slope)
        Ca = self.getCa(slope)
        Cl1 = self.getCl(slope)
        Cl2 = self.getCl(slope)
        Q1 = ((Ca+Cd)/2)**2
        Q2 = ((Cl1+Cl2)/2)**2
        D1 =  (Ca-Cd)/2
        D2 = (Cl2-Cl1)/2
        return getCAMIScost(beta,Q1,Q2,D1,D2)
    
    def getRawCost(self, slope, aspect, heading):
        beta = computeBeta(aspect,heading)
        Cd = self.getRawCd(slope)
        Ca = self.getRawCa(slope)
        Cl1 = self.getRawCl(slope)
        Cl2 = self.getRawCl(slope)
        Q1 = ((Ca+Cd)/2)**2
        Q2 = ((Cl1+Cl2)/2)**2
        D1 =  (Ca-Cd)/2
        D2 = (Cl2-Cl1)/2
        return getCAMIScost(beta,Q1,Q2,D1,D2)
    
    def showAnisotropy(self):
        linearGradient = np.linspace(0,89,90)
        heading = np.arange(0, 2*np.pi, 0.01)
        aspect = [1,0]
        Cd = np.zeros_like(linearGradient)
        Ca = np.zeros_like(linearGradient)
        Cl1 = np.zeros_like(linearGradient)
        Cl2 = np.zeros_like(linearGradient)
        for i,g in enumerate(linearGradient):
            Cd[i] = self.getCd(linearGradient[i])
            Ca[i] = self.getCa(linearGradient[i])
            Cl1[i] = self.getCl(linearGradient[i])
            Cl2[i] = self.getCl(linearGradient[i])
            
        Bs = []
        Cs = []
        Anisotropy = np.zeros_like(linearGradient)
        AnisotropyAD = np.zeros_like(linearGradient)
        AnisotropyAL = np.zeros_like(linearGradient)
        AnisotropyDL = np.zeros_like(linearGradient)
        for i,g in enumerate(linearGradient):
            for theta in heading:
                B = computeBeta(aspect,theta)
                Bs.append(np.arctan2(B[1],B[0]))
                preCost = computeCAMIScost(B,Cd[i],Ca[i],Cl1[i],Cl2[i])
                Cs.append(preCost)
            Anisotropy[i] = max(Cs)/min(Cs)
            AnisotropyAD[i] = Ca[i]/Cd[i]
            AnisotropyAL[i] = Ca[i]/Cl1[i]
            AnisotropyDL[i] = Cd[i]/Cl1[i]
            Cs = []
            
        fig, axes = plt.subplots()
        l1 = axes.plot(linearGradient,Anisotropy,color='b', linestyle='dashed', label = '$Anisotropy$')
        l2 = axes.plot(linearGradient,AnisotropyAD,color='r', linestyle='dashed', label = '$Ascent-Descent Anisotropy$')
        l3 = axes.plot(linearGradient,AnisotropyAL,color='g', linestyle='dashed', label = '$Ascent-Lateral Anisotropy$')
        l4 = axes.plot(linearGradient,AnisotropyDL,color='m', linestyle='dashed', label = '$Descent-Lateral Anisotropy$')
        lns = l1+l2+l3+l4
        labs = [l.get_label() for l in lns]
        axes.legend(lns, labs, fontsize='small')
        axes.set_xlim([0.0, 45.0])
        axes.set_ylim([0.0, 10.0])
        
    def showBraking(self):
        plt.style.use('seaborn-darkgrid')
        
        plt.rcParams["font.family"] = "Constantia"
        plt.rcParams['mathtext.fontset'] = 'cm'
        plt.rcParams['mathtext.rm'] = 'serif'
        steepnessArray = np.linspace(0.0,45.0,90+2)
        descentRR = np.zeros_like(steepnessArray)
        for i,steepness in enumerate(steepnessArray):
            descentRR[i] = self.getRRd(steepness)
        descentFunction = np.abs(self.friction - np.tan(steepnessArray*deg2rad))*self.kmg
        fig1, ax1 = plt.subplots(figsize=(5, 4),constrained_layout=True)
        ax1.plot(steepnessArray, descentFunction, linestyle = 'dashed', color = 'g')
        ax1.plot(steepnessArray, descentRR, color = 'g')
        ax1.plot(self.brakePoint01[0]*rad2deg, self.brakePoint01[1], 'o', color = 'g')
        ax1.plot(self.brakePoint02[0]*rad2deg, self.brakePoint02[1], 'o', color = 'g')
        ax1.plot(self.brakePoint[0]*rad2deg, self.brakePoint[1], 'o', color = 'g')
        ax1.set_xlabel('Steepness α [degrees]')
        ax1.set_ylabel('Energy [Ws/m]')
        ax1.legend(('|ρ - tan α|','$\mathbb{R}_d$'))
        ax1.annotate('$α_{b} = $' + '{0:.2f}'.format(self.brakePoint[0]*rad2deg) + ' degrees',
                    xy=(self.brakePoint[0]*rad2deg, 
                        self.brakePoint[1]),
                    xytext=(9, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='left', va='bottom')
        ax1.annotate('$α - Δ α_{b} = $' + '{0:.2f}'.format(self.brakePoint01[0]*rad2deg) + ' degrees',
                    xy=(self.brakePoint01[0]*rad2deg, 
                        self.brakePoint01[1]),
                    xytext=(9, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='left', va='bottom')
        ax1.annotate('$α + Δ α_{b} = $' + '{0:.2f}'.format(self.brakePoint02[0]*rad2deg) + ' degrees',
                    xy=(self.brakePoint02[0]*rad2deg, 
                        self.brakePoint02[1]),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='right', va='bottom')
        plt.style.use('default')
        
    def showDirCosts(self):
        plt.style.use('seaborn-darkgrid')
        
        plt.rcParams["font.family"] = "Constantia"
        plt.rcParams['mathtext.fontset'] = 'cm'
        plt.rcParams['mathtext.rm'] = 'serif'
        steepnessArray = np.linspace(0.0,45.0,90+2)*deg2rad
        
        ascentCAMIS = np.zeros_like(steepnessArray)
        descentCAMIS = np.zeros_like(steepnessArray)
        lateralCAMIS = np.zeros_like(steepnessArray)
        
        descentFunction = np.abs(self.friction - np.tan(steepnessArray))
        ascentFunction = self.friction + np.tan(steepnessArray)
        lateralFunction = self.friction*np.cos(steepnessArray)
        
        pitchRiskFunction = (steepnessArray - self.risk_pitch + self.risk_margin)/self.risk_margin*self.risk_gain
        rollRiskFunction = (steepnessArray - self.risk_roll + self.risk_margin)/self.risk_margin*self.risk_gain
        blockRiskFunction = (steepnessArray - self.risk_block + self.risk_margin)/self.risk_margin*self.risk_gain
        
        for i,steepness in enumerate(steepnessArray):
            ascentCAMIS[i] = self.getCa(steepness*rad2deg)
            descentCAMIS[i] = self.getCd(steepness*rad2deg)
            lateralCAMIS[i] = self.getCl(steepness*rad2deg)
        
        fig1, ax1 = plt.subplots(figsize=(5, 4),constrained_layout=True)
        ax1.plot(steepnessArray*rad2deg, ascentFunction, linestyle = 'dashed', color = 'b')
        ax1.plot(steepnessArray*rad2deg, descentFunction, linestyle = 'dashed', color = 'g')
        ax1.plot(steepnessArray*rad2deg, lateralFunction, linestyle = 'dashed', color = 'orange')
        ax1.plot(steepnessArray*rad2deg, ascentCAMIS, color = 'b')
        ax1.plot(steepnessArray*rad2deg, descentCAMIS, color = 'g')
        ax1.plot(steepnessArray*rad2deg, lateralCAMIS, color = 'orange')
#        ax1.plot(steepnessArray*rad2deg, blockRiskFunction, linestyle = 'dotted', color = 'y')
#        ax1.plot(steepnessArray*rad2deg, pitchRiskFunction, linestyle = 'dotted', color = 'r')
#        ax1.plot(steepnessArray*rad2deg, rollRiskFunction, linestyle = 'dotted', color = 'm')
#        ax1.plot(self.ascentBezierPoint_initial[0]*rad2deg,
#                 self.ascentBezierPoint_initial[1], 'o', color = 'b')
#        ax1.plot(self.descentBezierPoint_initial[0]*rad2deg,
#                 self.descentBezierPoint_initial[1], 'o', color = 'g')
#        ax1.plot(self.lateralBezierPoint_initial[0]*rad2deg,
#                 self.lateralBezierPoint_initial[1], 'o', color = 'orange')
#        ax1.plot(self.ascentBezierPoint_intersection1[0]*rad2deg,
#                 self.ascentBezierPoint_intersection1[1], 'o', color = 'b')
#        ax1.plot(self.ascentBezierPoint_intersection2[0]*rad2deg,
#                 self.ascentBezierPoint_intersection2[1], 'o', color = 'b')
#        ax1.plot(self.descentBezierPoint_intersection1[0]*rad2deg,
#                 self.descentBezierPoint_intersection1[1], 'o', color = 'g')
#        ax1.plot(self.descentBezierPoint_intersection2[0]*rad2deg,
#                 self.descentBezierPoint_intersection2[1], 'o', color = 'g')
#        ax1.plot(self.lateralBezierPoint_intersection1[0]*rad2deg,
#                 self.lateralBezierPoint_intersection1[1], 'o', color = 'orange')
#        ax1.plot(self.lateralBezierPoint_intersection2[0]*rad2deg,
#                 self.lateralBezierPoint_intersection2[1], 'o', color = 'orange')
#        ax1.plot(self.lateralBezierPoint_risk[0]*rad2deg,
#                 self.lateralBezierPoint_risk[1], 'o', color = 'orange')
#        ax1.plot(self.ascentBezierPoint_risk[0]*rad2deg,
#                 self.ascentBezierPoint_risk[1], 'o', color = 'b')
#        ax1.plot(self.brakePoint[0]*rad2deg, self.brakePoint[1], 'o', color = 'g')
#        ax1.plot(self.descentBezierPoint_risk[0]*rad2deg, 
#                 self.ascentBezierPoint_risk[1], 'o', color = 'g')
        ax1.legend(('ρ + tan α','|ρ - tan α|','ρ cos α','$R_a/ν_∥$','$R_d/ν_∥$','$R_l/ν_⟂$(1 + 1/ρ tan α)'))
        ax1.set_xlim([0.0, 35.0])
        ax1.set_ylim([0.0, self.getCa(35.0)])
        ax1.set_xlabel('Steepness α [degrees]')
        ax1.set_ylabel('Cost [Ws/m]')
        plt.style.use('default')
    def showCAMIS(self):
        
        linearGradient = np.linspace(0,89,90)
        heading = np.arange(0, 2*np.pi, 0.01)
        aspect = [1,0]
        
        fig = plt.figure()
        
        axes3 = fig.add_subplot(221, projection='3d')
        Xs = []
        Ys = []
        
        Cd = np.zeros_like(linearGradient)
        Ca = np.zeros_like(linearGradient)
        Cl1 = np.zeros_like(linearGradient)
        Cl2 = np.zeros_like(linearGradient)
         
        for i,g in enumerate(linearGradient):
            Cd[i] = self.getCd(linearGradient[i])
            Ca[i] = self.getCa(linearGradient[i])
            Cl1[i] = self.getCl(linearGradient[i])
            Cl2[i] = self.getCl(linearGradient[i])
            for theta in heading:
                B = computeBeta(aspect,theta)
                preCost = computeCAMIScost(B,Cd[i],Ca[i],Cl1[i],Cl2[i])
                Xs.append(B[0]*preCost)
                Ys.append(B[1]*preCost)
            axes3.plot(Xs, Ys, g, color=plt.cm.jet(float(g)/self.slopeThreshold))
            Xs = []
            Ys = []
        
        axes3.set_xlim(-self.Cmax,self.Cmax)
        axes3.set_ylim(-self.Cmax,self.Cmax)
        axes3.set_zlim(0,1.2*linearGradient[-1])
        axes3.set_xlabel('Parallel Cost',fontsize='medium')
        axes3.set_ylabel('Perpendicular Cost',fontsize='medium')
        axes3.set_zlabel('Steepness (deg)',fontsize='medium')
        axes3.xaxis.labelpad=-12
        axes3.yaxis.labelpad=-12
        axes3.zaxis.labelpad=-14
        axes3.view_init(elev=30, azim=-50)
        axes3.tick_params(axis="x",direction="in", pad=-6)
        axes3.tick_params(axis="y",direction="in", pad=-6)
        axes3.set_aspect('equal')
        axes3.dist = 7
        for spine in axes3.spines.values():
            spine.set_visible(False)
        
        
        axes5 = plt.subplot(222, projection='polar')
        axes6 = plt.subplot(224, projection='polar')
        axes5.set_facecolor('xkcd:light blue')
        axes6.set_facecolor('xkcd:light sage')
        Bs = []
        Ps = []
        Cs = []
        for i,g in enumerate(linearGradient):
            for theta in heading:
                B = computeBeta(aspect,theta)
                Bs.append(np.arctan2(B[1],B[0]))
                preCost = computeCAMIScost(B,Cd[i],Ca[i],Cl1[i],Cl2[i])
                Cs.append(preCost)
                Ps.append(1/preCost)
            axes5.plot(Bs, Cs, 'xkcd:sea blue', lw = 2, color = plt.cm.jet(float(g)/self.slopeThreshold))
            axes6.plot(Bs, Ps, 'xkcd:leaf', lw = 2, color = plt.cm.jet(float(g)/self.slopeThreshold))
            Cs = []
            Ps = []
            Bs = []
        fig.tight_layout()
        
        ax3 = fig.add_subplot(223)
        sqrtQ1 = np.zeros_like(linearGradient)
        sqrtQ2 = np.zeros_like(linearGradient)
        D1 = np.zeros_like(linearGradient)
        D2 = np.zeros_like(linearGradient)
        Cn = np.zeros_like(linearGradient)
        for i,g in enumerate(linearGradient):
            sqrtQ1[i] = (Ca[i]+Cd[i])/2
            sqrtQ2[i] = (Cl2[i]+Cl1[i])/2
            D1[i] = (Ca[i]-Cd[i])/2
            D2[i] = (Cl2[i]-Cl1[i])/2
#            Cn[i] = math.sqrt((Ca[i]+Cd[i])*(Cl1[i]+Cl2[i])/4)
        l1 = ax3.plot(linearGradient,sqrtQ1,color='b', label = '$K_∥$')
        l2 = ax3.plot(linearGradient,sqrtQ2,color='g', label = '$K_⊥$')
        l3 = ax3.plot(linearGradient,D1,color='m', label = '$D_∥$')
        #l4 = ax3.plot(linearGradient,D2,color='y', label = '$D_⊥$')
        l5 = ax3.plot(linearGradient,Ca,color='b', linestyle='dashed', label = '$C_a$')
        l6 = ax3.plot(linearGradient,Cd,color='g', linestyle='dashed', label = '$C_d$')
        l7 = ax3.plot(linearGradient,Cl1,color='m', linestyle='dashed', label = '$C_{l}$')
        #l8 = ax3.plot(linearGradient,Cl2,color='y', linestyle='dashed', label = '$C_{l2}$')
#        box = ax3.get_position()
#        ax3.set_position([box.x0, box.y0, box.width * 0.5, box.height])
        
        plt.ylabel('Cost [As/m]')
        plt.xlabel('Slope Gradient [degrees]')
        plt.grid('True')
        Bs = []
        Cs = []
        Anisotropy = np.zeros_like(linearGradient)
        for i,g in enumerate(linearGradient):
            for theta in heading:
                B = computeBeta(aspect,theta)
                Bs.append(np.arctan2(B[1],B[0]))
                preCost = computeCAMIScost(B,Cd[i],Ca[i],Cl1[i],Cl2[i])
                Cs.append(preCost)
            Anisotropy[i] = max(Cs)/min(Cs)
#            Cn[i] = math.sqrt((Ca[i]+Cd[i])*(Cl1[i]+Cl2[i])/4)/Cd[0]
            Cn[i] = .5*math.sqrt((Ca[i]*Cd[i]*(Cl1[i]**2+Cl2[i]**2)+Cl1[i]*Cl2[i]*(Ca[i]**2+Cd[i]**2))/math.sqrt(Cd[i]*Ca[i]*Cl1[i]*Cl2[i]))/Cd[0]
            Cs = []
            Bs = []
        ax3b = ax3.twinx()
        l9 = ax3b.plot(linearGradient,Anisotropy,color='r', label = '$ϒ$')
        l10 = ax3b.plot(linearGradient,Cn,color='r', linestyle='dashed', label = '$C_n/C_o$')
        ax3b.set_ylabel('Ratio', color='r')
        ax3b.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax3b.tick_params('y', colors='r')
        lns = l1+l2+l3+l5+l6 + l7 + l9 + l10
        labs = [l.get_label() for l in lns]
        ax3.legend(lns, labs, fontsize='small',\
            loc='upper right', bbox_to_anchor=(1.5, 1.1), ncol=1)
        fig.tight_layout()
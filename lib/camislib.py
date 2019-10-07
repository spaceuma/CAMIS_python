# -*- coding: utf-8 -*-

#============================CAMIS library=====================================
#           Continuous Anisotropic Model for Inclined Surfaces
#          Author: J. Ricardo Sanchez Ibanez (ricardosan@uma.es)
# -----------------------------------------------------------------------------
#   This file contains a library of python functions dedicated to the 
#==============================================================================

import numpy as np
from scipy.optimize import curve_fit
from bisect import bisect_left
import csv
import math

deg2rad = np.pi/180
rad2deg = 180/np.pi


# =============================================================================
#    rpy2ab -> Mapping from Roll-Pitch-Yaw angles to Alpha-Beta angles
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
        Gradient.append(rad2deg*np.arccos(np.cos(deg2rad*Roll[i])*np.cos(deg2rad*Pitch[i])))
        
    for i, r in enumerate(Roll):
        aspect = [np.cos(deg2rad*Yaw[i])*np.sin(deg2rad*Pitch[i])*np.cos(deg2rad*Roll[i]) + np.sin(deg2rad*Yaw[i])*np.sin(deg2rad*Roll[i]), np.sin(deg2rad*Yaw[i])*np.sin(deg2rad*Pitch[i])*np.cos(deg2rad*Roll[i]) - np.cos(deg2rad*Yaw[i])*np.sin(deg2rad*Roll[i])]
        aspect = aspect/np.linalg.norm(aspect)
        aspectVectorX.append(aspect[0])
        aspectVectorY.append(aspect[1])
        heading = [np.cos(deg2rad*Yaw[i])*np.cos(deg2rad*Pitch[i]), np.sin(deg2rad*Yaw[i])*np.cos(deg2rad*Pitch[i])]
        heading = heading/np.linalg.norm(heading)
        headingVectorX.append(heading[0])
        headingVectorY.append(heading[1])
        c, s = aspect[0], aspect[1]
        R = np.array(((c,s), (-s, c)))
        T = np.dot(R,np.transpose(heading))
        BetaX.append(T[0])
        BetaY.append(T[1])
        
    return Gradient, BetaX, BetaY

# =============================================================================
#    Explicit formulation of CAMIS to be used by a path planner
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
#    Bt = np.transpose(B)
#    Q = np.array(((Q1,D1*D2), (D1*D2, Q2)))
#    D = np.array(((D1), (D2)))
#    return (math.sqrt(np.dot(np.dot(Bt,Q),B))-np.dot(Bt,D))
    return np.sqrt(Q1*B[0]**2+Q2*B[1]**2+2*D1*D2*B[0]*B[1]) - \
           (B[0]*D1 + B[1]*D2)

# =============================================================================
#    Explicit formulation of CAMIS
# =============================================================================
           
           
    
def dirCost(gradient, K):
    # Gradient in degrees!
    return np.polyval(K, gradient)
    
def computeDirCosts(gradient, beta, cost):
    popt,_ = curve_fit(fittingCAMIS, (gradient,beta), cost)
    cdRoots = (popt[0],popt[1],popt[-1])
    caRoots = (8*popt[2],popt[3],popt[-1])
    cl1Roots = (popt[4],popt[5],popt[-1])
    cl2Roots = (popt[4],popt[5],popt[-1])
    return cdRoots, caRoots, cl1Roots, cl2Roots



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
    

def fittingCAMIS(x,x1,x2,x3,x4,x5, x6,Co):
    alpha, beta = x
    Cd = dirCost(alpha, [x1, x2, Co])
    Ca = dirCost(alpha, [x3, x4, Co])
    Cl1 = dirCost(alpha, [x5, x6, Co])
    Cl2 = dirCost(alpha, [x5, x6, Co])
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

def lookup(x, xs, ys):
    if x <= xs[0]:  return ys[0]
    if x >= xs[-1]: return ys[-1]

    i = bisect_left(xs, x)
    k = (x - xs[i-1])/(xs[i] - xs[i-1])
    y = k*(ys[i]-ys[i-1]) + ys[i-1]

    return y


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
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    linearGradient = np.linspace(0,30,31)
    heading = np.arange(0, 2*np.pi, 0.01)
    aspect = [1,0]
    
    Cd = dirCost(linearGradient, CdRoots)
    Ca = dirCost(linearGradient, CaRoots)
    Cl1 = dirCost(linearGradient, Cl1Roots)
    Cl2 = dirCost(linearGradient, Cl2Roots)
    
    CMax = max(Ca[-1],Cl1[-1])
    CMax = max(CMax,Cl2[-1])
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
    
    
    
    
    
    fig3 = plt.figure()
    
    sqrtQ1 = np.zeros_like(linearGradient)
    sqrtQ2 = np.zeros_like(linearGradient)
    D1 = np.zeros_like(linearGradient)
    D2 = np.zeros_like(linearGradient)
    for i,g in enumerate(linearGradient):
        sqrtQ1[i] = (Ca[i]+Cd[i])/2
        sqrtQ2[i] = (Cl2[i]+Cl1[i])/2
        D1[i] = (Ca[i]-Cd[i])/2
        D2[i] = (Cl2[i]-Cl1[i])/2
    ax3 = fig3.add_subplot(1, 1, 1)
    ax3.plot(linearGradient,sqrtQ1,color='b')
    ax3.plot(linearGradient,sqrtQ2,color='g')
    ax3.plot(linearGradient,D1,color='m')
    ax3.plot(linearGradient,D2,color='y')
    ax3.plot(linearGradient,Ca,color='burlywood')
    ax3.plot(linearGradient,Cd,color='lightblue')
    ax3.plot(linearGradient,Cl1,color='yellowgreen')
    ax3.plot(linearGradient,Cl2,color='coral')
    ax3.legend(['sqrt(Q1)','sqrt(Q2)','D1','D2','Ascent Cost Ca',\
                'Descent Cost Cd','Lateral Cost Cl1','Lateral Cost Cl2'])
    plt.ylabel('Cost [As/m]')
    plt.xlabel('Slope Gradient [degrees]')
    
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
        Cs = []
        Bs = []
    ax3b = ax3.twinx()
    ax3b.plot(linearGradient,Anisotropy,color='r')
    ax3b.set_ylabel('Anisotropy', color='r')
    ax3b.tick_params('y', colors='r')
    fig3.tight_layout()
    


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(Cd, np.zeros_like(linearGradient), linearGradient, color = 'r')
    ax.plot(-Ca, np.zeros_like(linearGradient), linearGradient, color = 'r')
    ax.plot(np.zeros_like(linearGradient), Cl1, linearGradient, color = 'r')
    ax.plot(np.zeros_like(linearGradient), -Cl2, linearGradient, color = 'r')
    Xs = []
    Ys = []
     
    for i,g in enumerate(linearGradient):
        for theta in heading:
            B = computeBeta(aspect,theta)
            preCost = computeCAMIScost(B,Cd[i],Ca[i],Cl1[i],Cl2[i])
            Xs.append(B[0]*preCost)
            Ys.append(B[1]*preCost)
        ax.plot(Xs, Ys, g, color=plt.cm.jet(float(g)/25))
        Xs = []
        Ys = []
    ax.set_xlim(-CMax,CMax)
    ax.set_ylim(-CMax,CMax)
    ax.set_aspect('equal')
    fig.tight_layout()
    
    fig5 = plt.figure()
    ax5 = fig5.add_subplot(1, 1, 1)
    ax5 = plt.subplot(111, projection='polar')
    ax5.set_facecolor('xkcd:light blue')
    Bs = []
    Cs = []
    for i,g in enumerate(linearGradient):
        for theta in heading:
            B = computeBeta(aspect,theta)
            Bs.append(np.arctan2(B[1],B[0]))
            preCost = computeCAMIScost(B,Cd[i],Ca[i],Cl1[i],Cl2[i])
            Cs.append(preCost)
        ax5.plot(Bs, Cs, 'xkcd:sea blue', lw = 2)
        Cs = []
        Bs = []
    fig5.tight_layout()
        
    fig6 = plt.figure()
    ax6 = fig6.add_subplot(1, 1, 1)
    ax6 = plt.subplot(111, projection='polar')
    ax6.set_facecolor('xkcd:light sage')
    Bs = []
    Ps = []
    for i,g in enumerate(linearGradient):
        for theta in heading:
            B = computeBeta(aspect,theta)
            Bs.append(np.arctan2(B[1],B[0]))
            preCost = computeCAMIScost(B,Cd[i],Ca[i],Cl1[i],Cl2[i])
            Ps.append(1/preCost)
        ax6.plot(Bs, Ps, 'xkcd:leaf', lw = 2)
        Ps = []
        Bs = []
    fig6.tight_layout()
    
    
    
    
    
    
    
        
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

class CamisModel:
    def __init__(self, cdRoots, caRoots, cl1Roots, cl2Roots, anicoLUT):
        self.cdRoots = cdRoots
        self.caRoots = caRoots
        self.cl1Roots = cl1Roots
        self.cl2Roots = cl2Roots
        self.anicoLUT = anicoLUT
        self.slopeThreshold = self.anicoLUT[0][-1]
    
    @classmethod
    def fromFile(cls, camis_file):
        return cls(*readCamis(camis_file))
    
    @classmethod
    def fromRoots(cls, cdRoots, caRoots, cl1Roots, cl2Roots, slopeThreshold):
        anicoLUT = computeAniCoLUT(cdRoots, caRoots, cl1Roots, cl2Roots,\
                                        slopeThreshold)
        return cls(cdRoots, caRoots, cl1Roots, cl2Roots, anicoLUT)
    
    def getCost(self, slope, aspect, heading):
        beta = computeBeta(aspect,heading)
        Cd = dirCost(slope, self.cdRoots)
        Ca = dirCost(slope, self.caRoots)
        Cl1 = dirCost(slope, self.cl1Roots)
        Cl2 = dirCost(slope, self.cl2Roots)
        Q1 = ((Ca+Cd)/2)**2
        Q2 = ((Cl1+Cl2)/2)**2
        D1 =  (Ca-Cd)/2
        D2 = (Cl2-Cl1)/2
        return getCAMIScost(beta,Q1,Q2,D1,D2)
        
    def printSlopeThreshold(self):
        print("The slope threshold is " + str(self.slopeThreshold) \
              + " degrees")
    def printCostRoots(self):
        print("Descent Cost Roots:   " + ''.join(str(c)+" " for c in \
                                               self.cdRoots))
        print("Ascent Cost Roots:    " + ''.join(str(c)+" " for c in \
                                               self.caRoots))
        print("Lateral 1 Cost Roots: " + ''.join(str(c)+" " for c in \
                                               self.cl1Roots))
        print("Lateral 2 Cost Roots: " + ''.join(str(c)+" " for c in \
                                               self.cl2Roots))
    def printAnicoLUT(self):
        print("Linear Gradient Array is:  " + ''.join(str(c)+" " for c in \
                                               self.anicoLUT[:][0]))
        print("Anisotropy Array is:       " + ''.join(str(c)+" " for c in \
                                               self.anicoLUT[:][1]))

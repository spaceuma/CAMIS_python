# -*- coding: utf-8 -*-

#============================CAMIS library=====================================
#           Continuous Anisotropic Model for Inclined Surfaces
#          Author: J. Ricardo Sanchez Ibanez (ricardosan@uma.es)
# -----------------------------------------------------------------------------
#   This file contains a library of python functions dedicated to the 
#==============================================================================

import numpy as np

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

# =============================================================================
#    Explicit formulation of CAMIS
# =============================================================================
def fittingCAMIS(x,x1,x2,x3,x4,x5,x6,x7,x8,Co):
    alpha, beta = x
    Cd = dirCost(alpha, x1, x2, Co)
    Ca = dirCost(alpha, x3, x4, Co)
    Cl1 = dirCost(alpha, x5, x6, Co)
    Cl2 = dirCost(alpha, x7, x8, Co)
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
    c, s = np.cos(aspect), np.sin(aspect)
    R = np.array(((c,s), (-s, c)))
    aHeading = np.array(((np.cos(heading)), (np.sin(heading))))
    T = np.dot(R,np.transpose(aHeading))
    return T

def dirCost(gradient, K1, K2, Co):
    # Gradient in degrees!
    return K1*gradient**2 + K2*gradient + Co

def getVectorialCostMap(slopeMap,popt):
    vectorialCostMap = np.zeros([4,slopeMap.shape[0],slopeMap.shape[1]])
    for i in range(slopeMap.shape[1]):
        for j in range(slopeMap.shape[0]):
            Cd = dirCost(slopeMap[j][i], popt[0], popt[1], popt[8])
            Ca = dirCost(slopeMap[j][i], popt[2], popt[3], popt[8])
            Cl1 = dirCost(slopeMap[j][i], popt[4], popt[5], popt[8])
            Cl2 = dirCost(slopeMap[j][i], popt[6], popt[7], popt[8])
            vectorialCostMap[0][j][i] = ((Ca+Cd)/2)**2 #Q1 = ((Ca+Cd)/2)**2
            vectorialCostMap[1][j][i] = ((Cl1+Cl2)/2)**2#Q2 = ((Cl2+Cl1)/2)**2
            vectorialCostMap[2][j][i] =  (Ca-Cd)/2 #D1 = (Ca-Cd)/2
            vectorialCostMap[3][j][i] = (Cl2-Cl1)/2 #D2 = (Cl2-Cl1)/2
    return vectorialCostMap

def showCAMIS(popt,beta,gradient,cost):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    linearGradient = np.linspace(0,30,31)
    heading = np.arange(0, 2*np.pi, 0.01)
    aspect = 0 
    
    Cd = dirCost(linearGradient, popt[0], popt[1], popt[8])
    Ca = dirCost(linearGradient, popt[2], popt[3], popt[8])
    Cl1 = dirCost(linearGradient, popt[4], popt[5], popt[8])
    Cl2 = dirCost(linearGradient, popt[6], popt[7], popt[8])
    
    CMax = max(Ca[-1],Cl1[-1])
    CMax = max(CMax,Cl2[-1])
    costX = np.multiply(np.cos(beta),cost)
    costY = np.multiply(np.sin(beta),cost)
    
    rgba_colors = np.zeros((costX.shape[0],4))
    rgba_colors[:,0] = np.abs(np.sin(beta))
    rgba_colors[:, 3] = np.abs(np.sin(beta))
    
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(1, 1, 1)
    ax1.set_facecolor('xkcd:pale grey')
    ax1.set_title('$\longleftarrow$ Lateral Direction #2      Lateral Direction #1 $\longrightarrow$')
    plt.grid('True')
    ax1.fill_betweenx(linearGradient,0,Cl1,color='xkcd:sea blue')
    ax1.fill_betweenx(linearGradient,0,-Cl2,color='xkcd:leaf')
    ax1.scatter(costY,gradient,color=rgba_colors)
    ax1.legend(['Lateral Cost 1','Lateral Cost 2','Experimental Data'])
    ax1.set_xlim(-CMax,CMax)
    ax1.set_ylim(0,linearGradient[-1])
    ax1.plot(Cl1,linearGradient, lw = 2, color = 'xkcd:dusk blue')
    ax1.plot(-Cl2,linearGradient, lw = 2, color = 'xkcd:camo green')
    plt.xlabel('Cost [Am/s]')
    plt.ylabel('Slope Gradient [degrees]')
    fig1.tight_layout()
    
    rgba_colors[:,0] = np.abs(np.cos(beta))
    rgba_colors[:, 3] = np.abs(np.cos(beta))
    
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1, 1, 1)
    ax2.set_facecolor('xkcd:pale grey')
    ax2.set_title('$\longleftarrow$ Ascent Direction      Descent Direction $\longrightarrow$')
    plt.grid('True')
    ax2.fill_betweenx(linearGradient,0,Cd,color='xkcd:sea blue')
    ax2.fill_betweenx(linearGradient,0,-Ca,color='xkcd:leaf')
    ax2.scatter(costX,gradient,color=rgba_colors)
    ax2.legend(['Descent Cost','Ascent Cost','Experimental Data'])
    ax2.set_xlim(-CMax,CMax)
    ax2.set_ylim(0,linearGradient[-1])
    ax2.plot(Cd,linearGradient, lw = 2, color = 'xkcd:dusk blue')
    ax2.plot(-Ca,linearGradient, lw = 2, color = 'xkcd:camo green')
    plt.xlabel('Cost [Am/s]')
    plt.ylabel('Slope Gradient [degrees]')
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
    ax3.legend(['sqrt(Q1)','sqrt(Q2)','D1','D2'])
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
        




# -*- coding: utf-8 -*-
"""
CAMIS
Continuous Anisotropic Model for Inclined Surfaces

This Python script 
"""

import numpy as np

def rpy2ab(Roll,Pitch,Yaw):
    deg2rad = np.pi/180
    rad2deg = 180/np.pi
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

def computeCAMIScost(heading,aspect,Cd,Ca,Cl1,Cl2):
    c, s = np.cos(aspect), np.sin(aspect)
    R = np.array(((c,s), (-s, c)))
    aHeading = np.array(((np.cos(heading)), (np.sin(heading))))
    T = np.dot(R,np.transpose(aHeading))
    Tt = np.transpose(T)
    K1 = (Ca+Cd)/2
    K2 = (Cl1+Cl2)/2
    D1 = (Ca-Cd)/2
    D2 = (Cl2-Cl1)/2
    Q1 = np.power(K1,2)
    Q2 = np.power(K2,2)
    Q3 = D1*D2
    Q = np.array(((Q1,Q3), (Q3, Q2)))
    D = np.array(((D1), (D2)))
    return (np.sqrt(np.dot(np.dot(Tt,Q),T))-np.dot(Tt,D)),T
        




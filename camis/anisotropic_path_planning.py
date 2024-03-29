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
"""

Python Implementation of bi-OUM

"""

import matplotlib.pyplot as plt
import numpy as np
import operator
import time
import bisect
import sys
from scipy.optimize import minimize_scalar
from numba import jit, float64, int64
from time import time
import math


# =============================================================================
#    Computation of the Total Cost Map
# =============================================================================
def computeTmap(VCmap,aspectMap,anisotropyMap,goal,start,Xmap,Ymap,res):
    # State Maps
    #  State -1 = far
    #  State 0 = narrow
    #  State 1 = accepted
    #  State 2 = closed
    stateMap = -1*np.ones_like(anisotropyMap)
    
    # Define nodeTargets
    nodeTarget = [goal[0],goal[1]]
    stateMap[nodeTarget[1],nodeTarget[0]] = 1
    
    # Both closedMaps are created
    stateMap = np.ones_like(anisotropyMap)*(-1)
    stateMap[nodeTarget[1],nodeTarget[0]] = 1
    stateMap[np.where(anisotropyMap == np.inf)] = 2
    stateMap[np.where(np.isnan(anisotropyMap))] = 2
    
    # Direction Maps are initialized
    dirMap = np.ones_like(anisotropyMap)*np.inf
    
    # Narrow Bands are initialized
    # Narrow Bands are initialized
    nbT = []
    nbNodes = []
    # Both Tmaps are initialized
    Tmap = np.ones_like(anisotropyMap)*np.inf
    Tmap[nodeTarget[1],nodeTarget[0]] = 0
    
    # Initial T update
    Tmap, nbT, nbNodes = updateNeighbours(nodeTarget, nbT, nbNodes, dirMap, Tmap, stateMap, VCmap, aspectMap, anisotropyMap, Xmap, Ymap, res, goal)

    while nbNodes:
        nodeTarget, nbT, nbNodes = getMinNB(nbT, nbNodes)
        stateMap[nodeTarget[1],nodeTarget[0]] = 1
        Tmap, nbT, nbNodes = updateNeighbours(nodeTarget, nbT, nbNodes, dirMap, Tmap, stateMap, VCmap, aspectMap, anisotropyMap, Xmap, Ymap, res)
        if stateMap[start[1],start[0]] == 2:
            break
    return Tmap, dirMap, stateMap


def computeBiTmap(VCmap,aspectMap,anisotropyMap,goal,start,Xmap,Ymap,res):
    VCmapS = np.ones_like(VCmap)
    VCmapS[0] = VCmap[0]
    VCmapS[1] = VCmap[1]
    VCmapS[2] = -VCmap[2]
    VCmapS[3] = -VCmap[3]
    
    # State Maps
    #  State -1 = far
    #  State 0 = narrow
    #  State 1 = accepted
    #  State 2 = closed
    stateMapG = -1*np.ones_like(anisotropyMap)
    stateMapS = stateMapG
    
    # Define nodeTargets
    nodeTargetG = [goal[0],goal[1]]
    nodeTargetS = [start[0],start[1]]
    stateMapG[nodeTargetG[1],nodeTargetG[0]] = 1
    stateMapS[nodeTargetS[1],nodeTargetS[0]] = 1
    
    # Both closedMaps are created
    stateMapG = np.ones_like(anisotropyMap)*(-1)
    stateMapG[nodeTargetG[1],nodeTargetG[0]] = 1
    stateMapG[np.where(anisotropyMap == np.inf)] = 2
    stateMapG[np.where(np.isnan(anisotropyMap))] = 2
    stateMapS = np.ones_like(anisotropyMap)*(-1)
    stateMapS[nodeTargetS[1],nodeTargetS[0]] = 1
    stateMapS[np.where(anisotropyMap == np.inf)] = 2
    stateMapS[np.where(np.isnan(anisotropyMap))] = 2
    
    # Direction Maps are initialized
    dirMapG = np.ones_like(anisotropyMap)*np.inf
    dirMapS = np.ones_like(anisotropyMap)*np.inf
    
    # Narrow Bands are initialized
    nbTG = []
    nbNodesG = []
    nbTS = []
    nbNodesS = []
    # Both Tmaps are initialized
    TmapG = np.ones_like(anisotropyMap)*np.inf
    TmapG[nodeTargetG[1],nodeTargetG[0]] = 0
    TmapS = np.ones_like(anisotropyMap)*np.inf
    TmapS[nodeTargetS[1],nodeTargetS[0]] = 0
    
    # Initial T update
    
    TmapG, nbTG, nbNodesG = updateNeighbours(nodeTargetG, nbTG, nbNodesG, dirMapG, TmapG, stateMapG, VCmap, aspectMap, anisotropyMap, Xmap,Ymap,res, goal)
    TmapS, nbTS, nbNodesS = updateNeighbours(nodeTargetS, nbTS, nbNodesS, dirMapS, TmapS, stateMapS, VCmapS, aspectMap, anisotropyMap, Xmap,Ymap,res, start)   
 
    nodeLink = []
    while nbNodesG or nbNodesS:
        if nbNodesG:
            nodeTargetG, nbTG, nbNodesG = getMinNB(nbTG, nbNodesG)
            stateMapG[nodeTargetG[1],nodeTargetG[0]] = 1
            TmapG, nbTG, nbNodesG = updateNeighbours(nodeTargetG, nbTG, nbNodesG, dirMapG, TmapG, stateMapG, VCmap, aspectMap, anisotropyMap, Xmap,Ymap,res)
       
        if nbNodesS:
            nodeTargetS, nbTS, nbNodesS = getMinNB(nbTS, nbNodesS)
            stateMapS[nodeTargetS[1],nodeTargetS[0]] = 1
            TmapS, nbTS, nbNodesS = updateNeighbours(nodeTargetS, nbTS, nbNodesS, dirMapS, TmapS, stateMapS, VCmapS, aspectMap, anisotropyMap, Xmap,Ymap,res)   

        if len(nodeLink) == 0:
            if stateMapS[nodeTargetG[1],nodeTargetG[0]] == 1:        
                d1 = dirMapS[nodeTargetG[1],nodeTargetG[0]]
                d2 = dirMapG[nodeTargetG[1],nodeTargetG[0]] + np.pi
                if np.arccos(np.cos(d1)*np.cos(d2)+np.sin(d1)*np.sin(d2)) < .2:
                    nodeLink = nodeTargetG
                    break
            if stateMapG[nodeTargetS[1],nodeTargetS[0]] == 1:
                nodeLink = nodeTargetS
                d1 = dirMapS[nodeTargetS[1],nodeTargetS[0]]
                d2 = dirMapG[nodeTargetS[1],nodeTargetS[0]] + np.pi
                if np.arccos(np.cos(d1)*np.cos(d2)+np.sin(d1)*np.sin(d2)) < .2:
                    nodeLink = nodeTargetS
                    break
        elif stateMapS[nodeLink[1],nodeLink[0]] == 2 and stateMapG[nodeLink[1],nodeLink[0]] == 2:
            break
#        if stateMapS[nodeTargetG[1],nodeTargetG[0]] == 2:

#        if stateMapG[nodeTargetS[1],nodeTargetS[0]] == 2:
#            d1 = dirMapS[nodeTargetS[1],nodeTargetS[0]]
#            d2 = dirMapG[nodeTargetS[1],nodeTargetS[0]] + np.pi
#            if np.arccos(np.cos(d1)*np.cos(d2)+np.sin(d1)*np.sin(d2)) < .2:
#                nodeLink = nodeTargetS
#                break
    return TmapG, TmapS, dirMapG, dirMapS, nodeLink, stateMapG, stateMapS, d1, d2 - np.pi

def updateNeighbours(nodeTarget, nbT, nbNodes, dirMap, Tmap, stateMap, VCmap, aspectMap, anisotropyMap, Xmap,Ymap,res, startingNode = []):
    NN = getNeighbours(nodeTarget)
    NN.append(nodeTarget)
    N = []
    NClist = []
    
    for i in range(len(NN)):
        # Detecting Inner Accepted Nodes
        if stateMap[NN[i][1],NN[i][0]] == 1:
            N = getNotAccepted(getNeighbours(NN[i]),stateMap)
            if len(N) == 0:
                stateMap[NN[i][1],NN[i][0]] = 2
        # Getting new Considered Nodes
        if stateMap[NN[i][1],NN[i][0]] == -1:
            stateMap[NN[i][1],NN[i][0]] = 0
            NClist.append(NN[i])
            
    for i in range(len(NClist)):
        nodeTarget = NClist[i]
        aspect = aspectMap[:,nodeTarget[1],nodeTarget[0]]
        anisotropy = anisotropyMap[nodeTarget[1],nodeTarget[0]]
        nfPairs = []
        if len(startingNode)!= 0:
            nfPairs = np.concatenate((startingNode,startingNode))
        else:
            if anisotropy == 1.0:
                R = int(2)
            else:
#                R = int(np.ceil(anisotropy*1.1547005383792517)*2 + 1)#2/sqrt(3)
                R = int(np.ceil(anisotropy))
            afList = []
#            for j in range(-R,R+1):
#                for k in range(-R,R+1):
#                    try:
#                        if stateMap[nodeTarget[1]+j,nodeTarget[0]+k]==1 and \
#                        computeDistance(nodeTarget+[k,j],nodeTarget,res) <= 2*res*R:
#                            afList.append([nodeTarget[0]+k,nodeTarget[1]+j])
#                    except:
#                        pass
#            R = 2
            for j in range(-R,1):
                for k in range(-R-j,R+1):
#                    print(k,j)
                    try:
                        if stateMap[nodeTarget[1]+j,nodeTarget[0]+k]==1:
                            afList.append([nodeTarget[0]+k,nodeTarget[1]+j])
                    except:
                        pass
            for j in range(1,R+1):
                for k in range(-R,R-j+1):
#                    print(k,j)
                    try:
                        if stateMap[nodeTarget[1]+j,nodeTarget[0]+k]==1:
                            afList.append([nodeTarget[0]+k,nodeTarget[1]+j])
                    except:
                        pass
            localAFPairs = []
            SS = []
            SS[:] = afList
            while (len(SS)!=0):
                ss = SS[0]
                del SS[0]
                for j in SS:
                    if (computeDistance(ss,j,res) <= 1.1*res):
                        localAFPairs.append(np.concatenate((ss,j)))
            nfPairs = []
            for j in localAFPairs:
                if checkNF(j,nodeTarget,anisotropy,res):
                    nfPairs.append(j)
            if len(nfPairs) == 0:
                for j in afList:
                    nfPairs.append(np.concatenate((j,j)))
        Q1 = VCmap[0,nodeTarget[1],nodeTarget[0]]
        Q2 = VCmap[1,nodeTarget[1],nodeTarget[0]]
        D1 = VCmap[2,nodeTarget[1],nodeTarget[0]]
        D2 = VCmap[3,nodeTarget[1],nodeTarget[0]]
        T,direc = computeT(nodeTarget, nfPairs, Q1, Q2, D1, D2, aspect, Tmap, dirMap,Xmap,Ymap,anisotropy,res)
        nIndex = bisect.bisect_left(nbT,T)
        nbT.insert(nIndex,T)
        nbNodes.insert(nIndex, nodeTarget)
        Tmap[nodeTarget[1],nodeTarget[0]] = T
        dirMap[nodeTarget[1],nodeTarget[0]] = direc
    return Tmap, nbT, nbNodes       

def getAccepted(nodeList, stateMap):
    acceptedList = []
    for i in range(len(nodeList)):  
        if stateMap[nodeList[i][1],nodeList[i][0]] == 1:
            acceptedList.append(nodeList[i])
    return acceptedList

def getNotAccepted(nodeList, stateMap):
    consideredList = []
    for i in range(len(nodeList)):     
        if stateMap[nodeList[i][1],nodeList[i][0]] <= 0:
            consideredList.append(nodeList[i])
    return consideredList

def getNeighbours(nodeTarget):
    nList = []
    neighbours = [[1,0],
                  [0,1],
                  [-1,1],
                  [-1,0],
                  [0,-1],
                  [1,-1],]
    for ni in neighbours:
        nN = np.add(nodeTarget,ni)
        nList.append(nN)
    return nList

def checkNF(afPair, n, anisotropy,res):
    C1 = afPair[2]-n[0];
    C2 = afPair[3]-n[1];
    C3 = afPair[0]-afPair[2];
    C4 = afPair[1]-afPair[3];
    amin = -.5*(2*C1*C3 + C1*C4 + C2*C3 + 2*C2*C4)
    ddp = 2*res**2
    dp0 = -ddp*amin
    p0 = .125*(3*ddp*C2**2 + ddp*(2*C1+C2)**2 -8*(res*anisotropy)**2)
    p1 = 0.5*ddp+dp0+p0
    dP = dp0**2 - 2*ddp*p0;
    if dP < 0:
        return False
    elif p0 <= 0 or p1 <= 0:
        return True
    elif amin >= 0 and amin <= 1:
        return True
    else:
        return False

def computeDistance(nodeA, nodeB, res):
    dx = res*(nodeA[0]+nodeA[1]/2 - nodeB[0] - nodeB[1]/2)
    dy = res*0.8660254037844386*(nodeA[1] - nodeB[1])
    return math.sqrt(dx**2+dy**2)

def computeT(nodeTarget, nfPairs, Q1, Q2, D1, D2, aspect, Tmap, dirMap,Xmap,Ymap,anisotropy,res):
    if np.isnan(aspect[0]) or np.isnan(aspect[1]):
        aspect = [1,0]
    T = np.inf
    direc = np.nan
    x = np.asarray([Xmap[nodeTarget[1],nodeTarget[0]],Ymap[nodeTarget[1],nodeTarget[0]]])
    if len(np.shape(nfPairs)) == 1:
        xj = np.asarray([Xmap[nfPairs[1],nfPairs[0]],Ymap[nfPairs[1],nfPairs[0]]])
        xk = np.asarray([Xmap[nfPairs[3],nfPairs[2]],Ymap[nfPairs[3],nfPairs[2]]])
        Tj = Tmap[nfPairs[1],nfPairs[0]]
        Tk = Tmap[nfPairs[3],nfPairs[2]]
        preT, preDir = optimizeCost(x,xj,xk,Tj,Tk,Q1,Q2,D1,D2,aspect,anisotropy)
        if preT < T:
            T = preT
            direc = preDir
    else:
        for i in range(len(nfPairs)):
            xj = np.asarray([Xmap[nfPairs[i][1],nfPairs[i][0]],Ymap[nfPairs[i][1],nfPairs[i][0]]])
            xk = np.asarray([Xmap[nfPairs[i][3],nfPairs[i][2]],Ymap[nfPairs[i][3],nfPairs[i][2]]])
            Tj = Tmap[nfPairs[i][1],nfPairs[i][0]]
            Tk = Tmap[nfPairs[i][3],nfPairs[i][2]]
            if ((anisotropy == 1.0)and(np.linalg.norm(x-xj) < 1.1*res)and(np.linalg.norm(x-xk) < 1.1*res)):
                preT, preDir = getEikonalCost(x,xj,xk,Tj,Tk,Q1)
            else:
                preT, preDir = optimizeCost(x,xj,xk,Tj,Tk,Q1,Q2,D1,D2,aspect,anisotropy)
            if preT < T:
                T = preT
                direc = preDir
    return T,direc
       
def optimizeCost(x,xj,xk,Tj,Tk,Q1,Q2,D1,D2,aspect,anisotropy):
    a = 0
    b = 1
    nIter = 500
    tau= 0.6180339887498949#(math.sqrt(5)-1)/2;
    epsilon1 = a + (1-tau)*(b-a)
    epsilon2 = a + tau*(b-a)
           
#    R = np.array(((aspect[0],aspect[1]), (-aspect[1], aspect[0])))
#    heading1 = x - epsilon1 * xj - (1-epsilon1) * xk
#    beta1 = np.dot(R,np.transpose(heading1)) 
    beta1x = aspect[0]*(x[0] - epsilon1 * xj[0] - (1-epsilon1) * xk[0]) + \
             aspect[1]*(x[1] - epsilon1 * xj[1] - (1-epsilon1) * xk[1])
    beta1y = -aspect[1]*(x[0] - epsilon1 * xj[0] - (1-epsilon1) * xk[0]) + \
              aspect[0]*(x[1] - epsilon1 * xj[1] - (1-epsilon1) * xk[1])   
    f_x1 = math.sqrt(Q1*beta1x**2+Q2*beta1y**2+2*D1*D2*beta1x*beta1y) + \
           (beta1x*D1 + beta1y*D2) + epsilon1*Tj + (1-epsilon1)*Tk
    
#    heading2 = x - epsilon2 * xj - (1-epsilon2) * xk
#    beta2 = np.dot(R,np.transpose(heading2))
    beta2x = aspect[0]*(x[0] - epsilon2 * xj[0] - (1-epsilon2) * xk[0]) + \
             aspect[1]*(x[1] - epsilon2 * xj[1] - (1-epsilon2) * xk[1])
    beta2y = -aspect[1]*(x[0] - epsilon2 * xj[0] - (1-epsilon2) * xk[0]) + \
              aspect[0]*(x[1] - epsilon2 * xj[1] - (1-epsilon2) * xk[1]) 
    f_x2 = math.sqrt(Q1*beta2x**2+Q2*beta2y**2+2*D1*D2*beta2x*beta2y) + \
           (beta2x*D1 + beta2y*D2) + epsilon2*Tj + (1-epsilon2)*Tk
    accuracy = 0.0001

    if (anisotropy > 10):
        nIter = 500/anisotropy
        accuracy = 0.0001/anisotropy
    k = 0

    while np.abs(b-a)>accuracy and k<nIter:
        k = k+1
        if f_x1 < f_x2:
            b = epsilon2
            epsilon2 = epsilon1
            epsilon1 = a + (1-tau)*(b-a)
        else:
            a = epsilon1
            epsilon1 = epsilon2
            epsilon2 = a + tau*(b-a)

#        heading1 = x - epsilon1 * xj - (1-epsilon1) * xk
#        beta1 = np.dot(R,np.transpose(heading1))     
#        f_x1 = math.sqrt(Q1*beta1[0]**2+Q2*beta1[1]**2+2*D1*D2*beta1[0]*beta1[1]) + \
#           (beta1[0]*D1 + beta1[1]*D2) + epsilon1*Tj + (1-epsilon1)*Tk
        beta1x = aspect[0]*(x[0] - epsilon1 * xj[0] - (1-epsilon1) * xk[0]) + \
             aspect[1]*(x[1] - epsilon1 * xj[1] - (1-epsilon1) * xk[1])
        beta1y = -aspect[1]*(x[0] - epsilon1 * xj[0] - (1-epsilon1) * xk[0]) + \
              aspect[0]*(x[1] - epsilon1 * xj[1] - (1-epsilon1) * xk[1])   
        f_x1 = math.sqrt(Q1*beta1x**2+Q2*beta1y**2+2*D1*D2*beta1x*beta1y) + \
           (beta1x*D1 + beta1y*D2) + epsilon1*Tj + (1-epsilon1)*Tk   
#        heading2 = x - epsilon2 * xj - (1-epsilon2) * xk
#        beta2 = np.dot(R,np.transpose(heading2))
#        f_x2 = math.sqrt(Q1*beta2[0]**2+Q2*beta2[1]**2+2*D1*D2*beta2[0]*beta2[1]) + \
#           (beta2[0]*D1 + beta2[1]*D2) + epsilon2*Tj + (1-epsilon2)*Tk
        beta2x = aspect[0]*(x[0] - epsilon2 * xj[0] - (1-epsilon2) * xk[0]) + \
             aspect[1]*(x[1] - epsilon2 * xj[1] - (1-epsilon2) * xk[1])
        beta2y = -aspect[1]*(x[0] - epsilon2 * xj[0] - (1-epsilon2) * xk[0]) + \
              aspect[0]*(x[1] - epsilon2 * xj[1] - (1-epsilon2) * xk[1]) 
        f_x2 = math.sqrt(Q1*beta2x**2+Q2*beta2y**2+2*D1*D2*beta2x*beta2y) + \
           (beta2x*D1 + beta2y*D2) + epsilon2*Tj + (1-epsilon2)*Tk
    
    if f_x1 < f_x2:
        minEpsilon = epsilon1
        minT = f_x1
    else:
        minEpsilon = epsilon2
        minT = f_x2
    
    minDirVector = (x-minEpsilon*xj-(1-minEpsilon)*xk)
    minDir = math.atan2(minDirVector[1],minDirVector[0])
    return minT,minDir

def getEikonalCost(x,xj,xk,Tj,Tk,Q):
    # We assume it is an irregular grid!!
    h = np.linalg.norm(x - xj)
    if (np.linalg.norm(xj-xk) < 0.01*h): #It is the same node
        T = Tj + h*math.sqrt(Q)
        cdir = x - xj
    else:
        # Be careful, here Q is cost**2
        T = (Tj + Tk + math.sqrt((Tj+Tk)**2 + 3*h**2*Q - 4*(Tj**2+Tk**2-Tj*Tk)))/2
        P = np.array([(x - xj)/h,
                      (x - xk)/h])
        cdir = (np.linalg.inv(P)).dot(np.array([(T - Tj)/h,
                                              (T - Tk)/h]))
    dirCharacteristic = math.atan2(cdir[1],cdir[0])
    return T,dirCharacteristic
    
def getMinNB(nbT,nbNodes):
    nodeTarget = nbNodes.pop(0)
    del nbT[0]
    return nodeTarget, nbT, nbNodes


def getPath(dirMap, IJ2XY, XY2IJ, initWaypoint, endWaypoint, Xmin, Ymin, res):
    path = []
    u = []
    path.append(initWaypoint)
    init = time()
    while time() - init < 10:
        waypoint = path[-1]
        uij = np.zeros_like(waypoint,int)
        if (uij[0] == np.inf)or(uij[1]==np.inf):
            raise ValueError('Waypoint with infinite Total Cost')
        try:
            uij[0] = int(np.round(interpolatePoint(waypoint-[Xmin,Ymin],XY2IJ[0])))
            uij[1] = int(np.round(interpolatePoint(waypoint-[Xmin,Ymin],XY2IJ[1])))
        except:
            print('ERROR: path is not fully computed')
            return path, u
        k1 = .1*res*interpolatedControl(waypoint,dirMap,uij,IJ2XY,res)
        u.append(uij)
        if any(np.isnan(k1)):
            break
        if (k1[0] == 0)and(k1[1] == 0):
            break
        waypoint = path[-1]-k1
#        k2 = .1*res*interpolatedControl(waypoint,dirMap,uij,IJ2XY,res)
#        if (k2[0] == 0)and(k2[1] == 0):
#            break
#        if any(np.isnan(k1)) or any(np.isnan(k2)):
#            break
#        if np.arccos((k1[0]*k2[0]+k1[1]*k2[1])/(np.linalg.norm(k1)*np.linalg.norm(k2)))>45.0*np.pi/180.0:
##            waypoint = path[-1]-10*k1
#            break
#        else:
#            waypoint = path[-1]-.5*(k1+k2)
        path.append(waypoint)
        if math.sqrt((path[-1][0] - endWaypoint[0])**2+(path[-1][1] - endWaypoint[1])**2) < 1.5*res:
            break
#        if (np.abs((k1[0] + k2[0])) <= np.finfo(float).eps) and (np.abs((k1[1] + k2[1])) <= np.finfo(float).eps):
#            return path, u
#    path.append(endWaypoint)
    return path, u

def interpolatedControl(waypoint,dirMap,uij,IJ2XY,res):
    xc = IJ2XY[:,uij[1],uij[0]]
    u,v,w,vij,wij = findSimplex(waypoint,xc,uij,IJ2XY,res)
    e1,e2,e3 = interpolationCoefficients(waypoint,u,v,w,res)
    if np.isnan(dirMap[uij[1],uij[0]])or(np.isinf(dirMap[uij[1],uij[0]])):
        d1 = [0,0]
    else:
        d1 = [np.cos(dirMap[uij[1],uij[0]]), np.sin(dirMap[uij[1],uij[0]])]
    if np.isnan(dirMap[vij[1],vij[0]])or(np.isinf(dirMap[vij[1],vij[0]])):
        d2 = [0,0]
    else:
        d2 = [np.cos(dirMap[vij[1],vij[0]]), np.sin(dirMap[vij[1],vij[0]])]
    if (np.isnan(dirMap[wij[1],wij[0]]))or(np.isinf(dirMap[wij[1],wij[0]])):
        d3 = [0,0]
    else:
        d3 = [np.cos(dirMap[wij[1],wij[0]]), np.sin(dirMap[wij[1],wij[0]])]
    
    control = [0,0]
    control[0] = e1*d1[0] + e2*d2[0] + e3*d3[0]
    control[1] = e1*d1[1] + e2*d2[1] + e3*d3[1]
    
    if (np.isnan(control[0]))or(np.isnan(control[1])):
        print()
    
    return control/np.linalg.norm(control)

def getTriInterpolation(waypoint,Tmap,XY2IJ,IJ2XY,res, Xmin, Ymin):
    uij = np.zeros_like(waypoint,int)
    uij[0] = int(np.round(interpolatePoint(waypoint-[Xmin,Ymin],XY2IJ[0])))
    uij[1] = int(np.round(interpolatePoint(waypoint-[Xmin,Ymin],XY2IJ[1])))
    xc = IJ2XY[:,uij[1],uij[0]]
    u,v,w,vij,wij = findSimplex(waypoint,xc,uij,IJ2XY,res)
    e1,e2,e3 = interpolationCoefficients(waypoint,u,v,w,res)
    T1 = Tmap[uij[1],uij[0]]
    T2 = Tmap[vij[1],vij[0]]
    T3 = Tmap[wij[1],wij[0]]
    return e1*T1 + e2*T2 + e3*T3

def getTriNearest(waypoint,Tmap,XY2IJ,IJ2XY,res, Xmin, Ymin):
    uij = np.zeros_like(waypoint,int)
    uij[0] = int(np.round(interpolatePoint(waypoint-[Xmin,Ymin],XY2IJ[0])))
    uij[1] = int(np.round(interpolatePoint(waypoint-[Xmin,Ymin],XY2IJ[1])))
    xc = IJ2XY[:,uij[1],uij[0]]
    u,v,w,vij,wij = findSimplex(waypoint,xc,uij,IJ2XY,res)
    e1,e2,e3 = interpolationCoefficients(waypoint,u,v,w,res)
    if e1 > e2:
        if e1 > e3:
            return Tmap[uij[1],uij[0]]
        else:
            return Tmap[wij[1],wij[0]]   
    else:
        if e2 > e3:
            return Tmap[vij[1],vij[0]]
        else:
            return Tmap[wij[1],wij[0]]

def getTriHighest(waypoint,Tmap,XY2IJ,IJ2XY,res, Xmin, Ymin):
    uij = np.zeros_like(waypoint,int)
    uij[0] = int(np.round(interpolatePoint(waypoint-[Xmin,Ymin],XY2IJ[0])))
    uij[1] = int(np.round(interpolatePoint(waypoint-[Xmin,Ymin],XY2IJ[1])))
    xc = IJ2XY[:,uij[1],uij[0]]
    u,v,w,vij,wij = findSimplex(waypoint,xc,uij,IJ2XY,res)
    return np.amax((Tmap[uij[1],uij[0]],Tmap[vij[1],vij[0]],Tmap[wij[1],wij[0]]))

def getTriLowest(waypoint,Tmap,XY2IJ,IJ2XY,res, Xmin, Ymin):
    uij = np.zeros_like(waypoint,int)
    uij[0] = int(np.round(interpolatePoint(waypoint-[Xmin,Ymin],XY2IJ[0])))
    uij[1] = int(np.round(interpolatePoint(waypoint-[Xmin,Ymin],XY2IJ[1])))
    xc = IJ2XY[:,uij[1],uij[0]]
    u,v,w,vij,wij = findSimplex(waypoint,xc,uij,IJ2XY,res)
    return np.amin((Tmap[uij[1],uij[0]],Tmap[vij[1],vij[0]],Tmap[wij[1],wij[0]]))
 
def getTriMedian(waypoint,Tmap,XY2IJ,IJ2XY,res, Xmin, Ymin):
    uij = np.zeros_like(waypoint,int)
    uij[0] = int(np.round(interpolatePoint(waypoint-[Xmin,Ymin],XY2IJ[0])))
    uij[1] = int(np.round(interpolatePoint(waypoint-[Xmin,Ymin],XY2IJ[1])))
    xc = IJ2XY[:,uij[1],uij[0]]
    u,v,w,vij,wij = findSimplex(waypoint,xc,uij,IJ2XY,res)
    return np.median((Tmap[uij[1],uij[0]],Tmap[vij[1],vij[0]],Tmap[wij[1],wij[0]]))
       
def interpolationCoefficients(x,u,v,w,h):
    C = 2/(3*h**2)
    z = x-w
    C1 = C*np.dot(z,u-v)
    C2 = C*np.dot(z,u-w)
    C3 = C*np.dot(z,v-w)
    e1 = C1 + C2
    e2 = C3 - C1
    e3 = 1 - C2 - C3
    return e1,e2,e3

def findSimplex(waypoint,xc,uij,IJ2XY,res):
    dVector = waypoint - xc
    ori = np.arctan2(dVector[1],dVector[0])*180/np.pi
    u = xc
    v = np.zeros_like(u)
    w = np.zeros_like(u)
    vij = np.zeros_like(uij,int)
    wij = np.zeros_like(uij,int)
    if ori > 120:
        vij[:] = uij + [-1,0]
        wij[:] = uij + [-1,1]
    elif ori > 60:
        vij[:] = uij + [-1,1]
        wij[:] = uij + [0,1]
    elif ori > 0:
        vij[:] = uij + [0,1]
        wij[:] = uij + [1,0]
    elif ori > -60:
        vij[:] = uij + [1,0]
        wij[:] = uij + [1,-1]
    elif ori > -120:
        vij[:] = uij + [1,-1]
        wij[:] = uij + [0,-1]
    else:
        vij[:] = uij + [0,-1]
        wij[:] = uij + [-1,0]
    v[:] = IJ2XY[:,vij[1],vij[0]]
    w[:] = IJ2XY[:,wij[1],wij[0]]
    return u,v,w,vij,wij
        
def interpolatePoint(point,mapI):
    i = np.uint32(np.fix(point[0]))
    j = np.uint32(np.fix(point[1]))
    a = point[0] - i
    b = point[1] - j
    
    
    m,n = np.uint32(mapI.shape)
    
    if i == n:
        if j == m:
            I = mapI[j,i]
        else:
            I = b*mapI[j+1,i] + (1-b)*mapI[j,i]
    else:
        if j == m:
            I = a*mapI[j,i+1] + (1-a)*mapI[j,i]
        else:
            a00 = mapI[j,i]
            a10 = mapI[j,i+1] - mapI[j,i]
            a01 = mapI[j+1,i] - mapI[j,i]
            a11 = mapI[j+1,i+1] + mapI[j,i] - mapI[j,i+1] - mapI[j+1,i]
            if a == 0:
                if b == 0:
                    I = a00
                else:
                    I = a00 + a01*b
            else:
                if b == 0:
                    I = a00 + a10*a
                else:
                    I = a00 + a10*a + a01*b + a11*a*b
                    
    return I

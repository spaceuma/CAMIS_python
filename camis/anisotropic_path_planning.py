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
def computeTmap(VCmap,aspectMap,anisotropyMap,goal,start,Xmap,Ymap,res,gridtype):
    # State Maps
    #  State -1 = far
    #  State 0 = narrow
    #  State 1 = accepted
    #  State 2 = closed
    stateMap = -1*np.ones_like(anisotropyMap)
    
    maxAnisoMap = np.ones_like(anisotropyMap)
    
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
    Tmap, nbT, nbNodes = updateNeighbours(nodeTarget, nbT, nbNodes, dirMap, Tmap, stateMap, VCmap, aspectMap, anisotropyMap, maxAnisoMap,  Xmap, Ymap, res, gridtype,goal)

    while nbNodes:
        nodeTarget, nbT, nbNodes = getMinNB(nbT, nbNodes)
        stateMap[nodeTarget[1],nodeTarget[0]] = 1
        Tmap, nbT, nbNodes = updateNeighbours(nodeTarget, nbT, nbNodes, dirMap, Tmap, stateMap, VCmap, aspectMap, anisotropyMap, maxAnisoMap,  Xmap, Ymap, res, gridtype)
        if stateMap[start[1],start[0]] == 2:
            break
    return Tmap, dirMap, stateMap


def computeBiTmap(VCmap, aspectMap, anisotropyMap, goal, start, Xmap, Ymap, 
                  res, gridtype):
    VCmapS = np.ones_like(VCmap)
    VCmapS[0] = VCmap[0]
    VCmapS[1] = VCmap[1]
    VCmapS[2] = -VCmap[2]
    VCmapS[3] = -VCmap[3]
    
    #This indicates the max anisotropy of any update the node is involved with
    maxAnisoMapG = np.ones_like(anisotropyMap)
    maxAnisoMapS = np.ones_like(anisotropyMap)
    
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
    
    TmapG, dirMapG, maxAnisoMapG, nbTG, nbNodesG = updateNeighbours(
        nodeTargetG, nbTG, nbNodesG, dirMapG, TmapG, stateMapG, VCmap, 
        aspectMap, anisotropyMap, maxAnisoMapG,  Xmap, Ymap, res, gridtype, 
        goal)
    TmapS, dirMapS, maxAnisoMapS, nbTS, nbNodesS = updateNeighbours(
        nodeTargetS, nbTS, nbNodesS, dirMapS, TmapS, stateMapS, VCmapS, 
        aspectMap, anisotropyMap, maxAnisoMapS, Xmap, Ymap, res, gridtype,
        start)   
 
    nodeLink = []
    while nbNodesG or nbNodesS:
        if nbNodesG:
            nodeTargetG, nbTG, nbNodesG = getMinNB(nbTG, nbNodesG)
            stateMapG[nodeTargetG[1],nodeTargetG[0]] = 1
            TmapG, dirMapG, maxAnisoMapG, nbTG, nbNodesG = updateNeighbours(
                nodeTargetG, nbTG, nbNodesG, dirMapG, TmapG, stateMapG, VCmap,  
                aspectMap, anisotropyMap, maxAnisoMapG, Xmap, Ymap, res, 
                gridtype)
            TmapG, dirMapG, maxAnisoMapG, nbTG, nbNodesG = updateTNarrowBand(
                nodeTargetG, nbTG, nbNodesG, dirMapG, TmapG, stateMapG, VCmap, 
                aspectMap, anisotropyMap, maxAnisoMapG, Xmap, Ymap, res, 
                gridtype)
       
        if nbNodesS:
            nodeTargetS, nbTS, nbNodesS = getMinNB(nbTS, nbNodesS)
            stateMapS[nodeTargetS[1],nodeTargetS[0]] = 1
            TmapS, dirMapS, maxAnisoMapS, nbTS, nbNodesS = updateNeighbours(
                nodeTargetS, nbTS, nbNodesS, dirMapS, TmapS, stateMapS, VCmapS, 
                aspectMap, anisotropyMap, maxAnisoMapS, Xmap, Ymap, res, 
                gridtype)  
            TmapS, dirMapS, maxAnisoMapS, nbTS, nbNodesS = updateTNarrowBand(
                nodeTargetS, nbTS, nbNodesS, dirMapS, TmapS, stateMapS, VCmapS, 
                aspectMap, anisotropyMap, maxAnisoMapS, Xmap, Ymap, res, 
                gridtype)

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
        elif stateMapS[nodeLink[1],nodeLink[0]] == 2 and stateMapG[nodeLink[1],
                                                             nodeLink[0]] == 2:
            break
    return TmapG, TmapS, dirMapG, dirMapS, nodeLink, stateMapG, stateMapS, \
           d1, d2 - np.pi




def updateNeighbours(nodeTarget, nbT, nbNodes, dirMap, Tmap, stateMap, VCmap, 
                     aspectMap, anisotropyMap, maxAnisoMap, Xmap, Ymap, res, 
                     gridtype, startingNode = []):
    NN = getNeighbours(nodeTarget, gridtype)
    NN.append(nodeTarget)
    N = []
    NClist = []
    
    for i in range(len(NN)):
        # Detecting Inner Accepted Nodes
        if stateMap[NN[i][1],NN[i][0]] == 1:
            N = getNotAccepted(getNeighbours(NN[i], gridtype), stateMap)
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
                R = int(1)
            else:
                R = int(np.ceil(anisotropy*1.1547005383792517))#2/sqrt(3)
            afList = []
            
            if gridtype == 'hex':
                afList = getHexAFlist(nodeTarget, R, stateMap, anisotropy, 
                                      afList)
            elif gridtype == 'sq':
                afList = getSqAFlist(nodeTarget, R, stateMap, anisotropy, 
                                      afList)
            else:
                raise ValueError('Wrong grid type')
                
            localAFPairs = []
            SS = []
            SS[:] = afList
            while (len(SS)!=0):
                ss = SS[0]
                del SS[0]
                if gridtype == 'sq':
                    localAFPairs.append(np.concatenate((ss,ss)))
                for j in SS:
                    if gridtype == 'hex' and areHexNeighbours(ss,j):
                        localAFPairs.append(np.concatenate((ss,j)))
                    if gridtype == 'sq' and areSq8Neighbours(ss,j):
                        localAFPairs.append(np.concatenate((ss,j)))
                        
            nfPairs = []
            for j in localAFPairs:
                if gridtype == 'hex' and checkHexNF(j,nodeTarget,anisotropy,
                                                    res):
                    nfPairs.append(j)
                if gridtype == 'sq' and checkSqNF(j,nodeTarget,anisotropy,
                                                    res):
                    nfPairs.append(j)
            
            if len(nfPairs) == 0:
                for j in afList:
                    nfPairs.append(np.concatenate((j,j)))
                    
        Q1 = VCmap[0,nodeTarget[1],nodeTarget[0]]
        Q2 = VCmap[1,nodeTarget[1],nodeTarget[0]]
        D1 = VCmap[2,nodeTarget[1],nodeTarget[0]]
        D2 = VCmap[3,nodeTarget[1],nodeTarget[0]]
        T,direc,nfAnisotropy = computeT(nodeTarget, nfPairs, Q1, Q2, D1, D2, 
                                        aspect, Tmap, dirMap, Xmap, Ymap,
                                        anisotropy, anisotropyMap, res,
                                        gridtype)
        nIndex = bisect.bisect_left(nbT,T)
        nbT.insert(nIndex,T)
        nbNodes.insert(nIndex, nodeTarget)
        Tmap[nodeTarget[1],nodeTarget[0]] = T
        dirMap[nodeTarget[1],nodeTarget[0]] = direc
        maxAnisoMap[nodeTarget[1],nodeTarget[0]] = nfAnisotropy
    return Tmap, dirMap, maxAnisoMap, nbT, nbNodes 

def getSqAFlist(nodeTarget, R, stateMap, anisotropy, afList):
    for j in range(-R,R+1):
        for k in range(-R,R+1):
            try:
                if stateMap[nodeTarget[1]+j,nodeTarget[0]+k]==1:
                    if anisotropy == 1.0:
                        if areSq8Neighbours(nodeTarget+[k,j],nodeTarget):
                            afList.append([nodeTarget[0]+k,nodeTarget[1]+j])
                    else:
                        afList.append([nodeTarget[0]+k,nodeTarget[1]+j])
            except:
                pass
    return afList


def getHexAFlist(nodeTarget, R, stateMap, anisotropy, afList):
    for j in range(-R,1):
        for k in range(-R-j,R+1):
            try:
                if stateMap[nodeTarget[1]+j,nodeTarget[0]+k]==1:
                    if anisotropy == 1.0:
                        if areHexNeighbours(nodeTarget+[k,j],nodeTarget):
                            afList.append([nodeTarget[0]+k,nodeTarget[1]+j])
                    else:
                        afList.append([nodeTarget[0]+k,nodeTarget[1]+j])
            except:
                pass
    for j in range(1,R+1):
        for k in range(-R,R-j+1):
            try:
                if stateMap[nodeTarget[1]+j,nodeTarget[0]+k]==1:
                    if anisotropy == 1.0:
                        if areHexNeighbours(nodeTarget+[k,j],nodeTarget):
                            afList.append([nodeTarget[0]+k,nodeTarget[1]+j])
                    else:
                        afList.append([nodeTarget[0]+k,nodeTarget[1]+j])
            except:
                pass 
    return afList


def getSqConsideredList(node, R, stateMap, relAnisotropy, consideredList):
    for j in range(-R,R+1):
        for k in range(-R,R+1):
            try:
                if stateMap[node[1]+j,node[0]+k]==0 and not any(
                        node[0]+k==x[0] and node[1]+j==x[1] for x in 
                        consideredList):
                    if relAnisotropy == 1.0:
                        if areSq8Neighbours(node+[k,j],node):
                            consideredList.append([node[0]+k,node[1]+j])
                    else:
                        consideredList.append([node[0]+k,node[1]+j])
            except:
                pass
    return consideredList

def getHexConsideredList(node, R, stateMap, relAnisotropy, consideredList):
    for j in range(-R,1):
        for k in range(-R-j,R+1):
            try:
                if stateMap[node[1]+j,node[0]+k]==0 and not any(
                        node[0]+k==x[0] and node[1]+j==x[1] for x in 
                        consideredList):
                    if relAnisotropy == 1.0:
                        if areHexNeighbours(node+[k,j],node):
                            consideredList.append([node[0]+k,node[1]+j])
                    else:
                        consideredList.append([node[0]+k,node[1]+j])
            except:
                pass
    for j in range(1,R+1):
        for k in range(-R,R-j+1):
            try:
                if stateMap[node[1]+j,node[0]+k]==0 and not any(
                        node[0]+k==x[0] and node[1]+j==x[1] for x in 
                        consideredList):
                    if relAnisotropy == 1.0:
                        if areHexNeighbours(node+[k,j],node):
                            consideredList.append([node[0]+k,node[1]+j])
                    else:
                        consideredList.append([node[0]+k,node[1]+j])
            except:
                pass 
    return consideredList

# updateTNarrowBand
# This function re-evaluates certain considered nodes
def updateTNarrowBand(nodeTarget, nbT, nbNodes, dirMap, Tmap, stateMap, VCmap, 
                      aspectMap, anisotropyMap, maxAnisoMap, Xmap, Ymap, res,
                      gridtype):
    
    #Initialization of lists
    consideredList = [] 
    
    # First, a set with neigh(NN)+NN is created
    NN = getNeighbours(nodeTarget, gridtype)
    NN.append(nodeTarget)
    # Only the accepted are obtained
    subAFlist = getAccepted(NN,stateMap)
    
    #If nodeTarget is isotropic, this is more simple
    if maxAnisoMap[nodeTarget[1],nodeTarget[0]] < 1.01:
        node = np.empty([2],dtype=int)
        node[:] = nodeTarget[:]
        relAnisotropy = maxAnisoMap[node[1],node[0]]
        R = int(1)
        if gridtype == 'hex':
            consideredList = getHexConsideredList(node, R, stateMap, 
                                                  relAnisotropy, 
                                                  consideredList)
        elif gridtype == 'sq':
            consideredList = getSqConsideredList(node, R, stateMap, 
                                                 relAnisotropy, consideredList)
        else:
            raise ValueError('Wrong grid type')
        
    else:
        #Now the considered nodes to re-evaluate are selected   
        for i in range(len(subAFlist)):
            node = np.empty([2],dtype=int)
            node[:] = subAFlist[i][:]
            relAnisotropy = maxAnisoMap[node[1],node[0]]
            # R = int(np.ceil(relAnisotropy) + 1)
            if relAnisotropy == 1.0:
                R = int(1)
            else:
                R = int(np.ceil(relAnisotropy*1.1547005383792517))#2/sqrt(3)
            if gridtype == 'hex':
                consideredList = getHexConsideredList(node, R, stateMap, 
                                                      relAnisotropy, 
                                                      consideredList)
            elif gridtype == 'sq':
                consideredList = getSqConsideredList(node, R, stateMap, 
                                                     relAnisotropy, 
                                                     consideredList)
            else:
                raise ValueError('Wrong grid type')
    
    for i in range(len(consideredList)):
        nodeTarget = np.empty([2],dtype=int)
        nodeTarget[:] = consideredList[i][:]
        nfPairs = []
        aspect = aspectMap[:,nodeTarget[1],nodeTarget[0]]
        anisotropy = anisotropyMap[nodeTarget[1],nodeTarget[0]]
        if anisotropy == 1.0:
            R = int(1)
        else:
            R = int(np.ceil(anisotropy*1.1547005383792517))#2/sqrt(3)
        afList = []
        if gridtype == 'hex':
            afList = getHexAFlist(nodeTarget, R, stateMap, anisotropy, 
                                  afList)
        elif gridtype == 'sq':
            afList = getSqAFlist(nodeTarget, R, stateMap, anisotropy, 
                                  afList)
        else:
            raise ValueError('Wrong grid type')
            
        localAFPairs = []
        SS = []
        SS[:] = afList
        while (len(SS)!=0):
            ss = SS[0]
            del SS[0]
            if gridtype == 'sq':
                localAFPairs.append(np.concatenate((ss,ss)))
            for j in SS:
                if gridtype == 'hex' and areHexNeighbours(ss,j):
                    localAFPairs.append(np.concatenate((ss,j)))
                if gridtype == 'sq' and areSq8Neighbours(ss,j):
                    localAFPairs.append(np.concatenate((ss,j)))
                    
        nfPairs = []
        for j in localAFPairs:
            if gridtype == 'hex' and checkHexNF(j,nodeTarget,anisotropy,
                                                res):
                nfPairs.append(j)
            if gridtype == 'sq' and checkSqNF(j,nodeTarget,anisotropy,
                                                res):
                nfPairs.append(j)
                
        if not len(nfPairs) == 0:
            Q1 = VCmap[0,nodeTarget[1],nodeTarget[0]]
            Q2 = VCmap[1,nodeTarget[1],nodeTarget[0]]
            D1 = VCmap[2,nodeTarget[1],nodeTarget[0]]
            D2 = VCmap[3,nodeTarget[1],nodeTarget[0]]
            T,direc,nfAnisotropy = computeT(nodeTarget, nfPairs, Q1, Q2, D1, D2, aspect, 
                               Tmap, dirMap,Xmap,Ymap,anisotropy,anisotropyMap,res,gridtype)
            if T < Tmap[nodeTarget[1],nodeTarget[0]]-0.0001:
                tempT = Tmap[nodeTarget[1],nodeTarget[0]]
                nIndex = bisect.bisect_left(nbT,tempT)
                try:
                    nIndex = next(x for x,n in enumerate(nbNodes[nIndex-1:],nIndex) if np.array_equal(nodeTarget,nbNodes[x]))
                    del nbT[nIndex]
                    del nbNodes[nIndex]
                    nIndex = bisect.bisect_left(nbT,T)
                    nbT.insert(nIndex,T)
                    nbNodes.insert(nIndex, nodeTarget)
                    Tmap[nodeTarget[1],nodeTarget[0]] = T  
                    dirMap[nodeTarget[1],nodeTarget[0]] = direc
                    maxAnisoMap[nodeTarget[1],nodeTarget[0]] = nfAnisotropy
                except IndexError:
                    print('What the...')
                    exit()
                except StopIteration:
                    print('What the...')
                    exit()
    return Tmap, dirMap, maxAnisoMap, nbT, nbNodes



      

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

def getNeighbours(nodeTarget, gridtype):
    nList = []
    if gridtype == 'hex':
        neighbours = [[1,0],
                      [0,1],
                      [-1,1],
                      [-1,0],
                      [0,-1],
                      [1,-1],]
    elif gridtype == 'sq':
        neighbours = [[1,0],
                      [0,1],
                      [-1,0],
                      [0,-1],]
    for ni in neighbours:
        nN = np.add(nodeTarget,ni)
        nList.append(nN)
    return nList

# This function returns whether two nodes are hexagonal neighbours or not
#@jit(nopython=True)
def areHexNeighbours(n1,n2):
    neighbours = [[1,0], [0,1], [-1,1], [-1,0], [0,-1], [1,-1],]
    dx = n1[0] - n2[0]
    dy = n1[1] - n2[1]
    return any(dx==x[0] and dy==x[1] for x in neighbours)

def areSqNeighbours(n1,n2):
    neighbours = [[1,0], [0,1], [-1,0], [0,-1],]
    dx = n1[0] - n2[0]
    dy = n1[1] - n2[1]
    return any(dx==x[0] and dy==x[1] for x in neighbours)

def areSq8Neighbours(n1,n2):
    neighbours = [[1,0], [0,1], [-1,0], [0,-1],
                  [-1,-1], [-1,1], [1,1], [1,-1],]
    dx = n1[0] - n2[0]
    dy = n1[1] - n2[1]
    return any(dx==x[0] and dy==x[1] for x in neighbours)

def checkSqNF(afPair, n, anisotropy,res):
    C1 = afPair[2]-n[0];
    C2 = afPair[3]-n[1];
    C3 = afPair[0]-afPair[2];
    C4 = afPair[1]-afPair[3];
    C5 = afPair[0]-n[0];
    C6 = afPair[1]-n[1];
    if anisotropy < 1.01:
        return math.sqrt(C1**2+C2**2) <= anisotropy and math.sqrt(
            C5**2+C6**2) <= anisotropy
        # return areSqNeighbours(afPair[0:2], n) and areSqNeighbours(
        #         afPair[2:4], n)
    else:
        C1 = afPair[2]-n[0];
        C2 = afPair[3]-n[1];
        C3 = afPair[0]-n[0];
        C4 = afPair[1]-n[1];
        return math.sqrt(C1**2+C2**2) <= anisotropy or math.sqrt(
            C3**2+C4**2) <= anisotropy



#@jit(nopython=True)  
def checkHexNF(afPair, n, anisotropy,res):
    C1 = afPair[2]-n[0];
    C2 = afPair[3]-n[1];
    C3 = afPair[0]-afPair[2];
    C4 = afPair[1]-afPair[3];
    # C5 = afPair[0]-n[0];
    # C6 = afPair[1]-n[1];
    if anisotropy < 1.01:
        return areHexNeighbours(afPair[0:2], n) and areHexNeighbours(
                afPair[2:4], n)
        # neighbours = [[1,0], [0,1], [-1,1], [-1,0], [0,-1], [1,-1],]
        # if any(C1==x[0] and C2==x[1] for x in neighbours) and any(
        #         C5==x[0] and C6==x[1] for x in neighbours):
        #     return True
        # else:
        #     return False
    else:
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

# @jit(nopython=True)  
# def computeDistance(nodeA, nodeB, res):
#     dx = res*(nodeA[0]+nodeA[1]/2 - nodeB[0] - nodeB[1]/2)
#     dy = res*0.8660254037844386*(nodeA[1] - nodeB[1])
#     return math.sqrt(dx**2+dy**2)

def computeT(nodeTarget, nfPairs, Q1, Q2, D1, D2, aspect, Tmap, dirMap, Xmap,
             Ymap, anisotropy, anisotropyMap, res, gridtype):
    if np.isnan(aspect[0]) or np.isnan(aspect[1]):
        aspect = [1,0]
    T = np.inf
    direc = np.nan
    nfAnisotropy = 1.0
    x = np.asarray([Xmap[nodeTarget[1],nodeTarget[0]],
                    Ymap[nodeTarget[1],nodeTarget[0]]])
    if len(np.shape(nfPairs)) == 1:
        xj = np.asarray([Xmap[nfPairs[1],nfPairs[0]],Ymap[nfPairs[1],nfPairs[0]]])
        xk = np.asarray([Xmap[nfPairs[3],nfPairs[2]],Ymap[nfPairs[3],nfPairs[2]]])
        Tj = Tmap[nfPairs[1],nfPairs[0]]
        Tk = Tmap[nfPairs[3],nfPairs[2]]
        # ToDo: create getEikonal for sq option
        if (anisotropy < 1.0001) and gridtype == 'hex':
            preT, preDir = getEikonalCost(x,xj,xk,Tj,Tk,Q1)
            # preT, preDir = optimizeCost(x,xj,xk,Tj,Tk,Q1,Q2,D1,D2,aspect,
            #                             anisotropy)
        else:
            preT, preDir = optimizeCost(x,xj,xk,Tj,Tk,Q1,Q2,D1,D2,aspect,
                                        anisotropy)
        if preT < T:
            T = preT
            direc = preDir
            nfAnisotropy = max(anisotropyMap[nfPairs[1],nfPairs[0]],
                               anisotropyMap[nfPairs[3],nfPairs[2]])
    else:
        for i in range(len(nfPairs)):
            xj = np.asarray([Xmap[nfPairs[i][1],nfPairs[i][0]],
                             Ymap[nfPairs[i][1],nfPairs[i][0]]])
            xk = np.asarray([Xmap[nfPairs[i][3],nfPairs[i][2]],
                             Ymap[nfPairs[i][3],nfPairs[i][2]]])
            Tj = Tmap[nfPairs[i][1],nfPairs[i][0]]
            Tk = Tmap[nfPairs[i][3],nfPairs[i][2]]
            # ToDo: create getEikonal for sq option
            if (anisotropy < 1.0001) and gridtype == 'hex':
                preT, preDir = getEikonalCost(x,xj,xk,Tj,Tk,Q1)
                # preT, preDir = optimizeCost(x,xj,xk,Tj,Tk,Q1,Q2,D1,D2,aspect,
                #                             anisotropy)
            else:
                preT, preDir = optimizeCost(x,xj,xk,Tj,Tk,Q1,Q2,D1,D2,aspect,
                                            anisotropy)
            if preT < T:
                T = preT
                direc = preDir
                nfAnisotropy = max(anisotropyMap[nfPairs[i][1],nfPairs[i][0]],
                                   anisotropyMap[nfPairs[i][3],nfPairs[i][2]])
    return T,direc, nfAnisotropy
   
@jit(nopython=True)  
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

        beta1x = aspect[0]*(x[0] - epsilon1 * xj[0] - (1-epsilon1) * xk[0]) + \
             aspect[1]*(x[1] - epsilon1 * xj[1] - (1-epsilon1) * xk[1])
        beta1y = -aspect[1]*(x[0] - epsilon1 * xj[0] - (1-epsilon1) * xk[0]) + \
              aspect[0]*(x[1] - epsilon1 * xj[1] - (1-epsilon1) * xk[1])   
        f_x1 = math.sqrt(Q1*beta1x**2+Q2*beta1y**2+2*D1*D2*beta1x*beta1y) + \
           (beta1x*D1 + beta1y*D2) + epsilon1*Tj + (1-epsilon1)*Tk   
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

# @jit(nopython=True)  
def getEikonalCost(x,xj,xk,Tj,Tk,Q):
    # We assume it is a regular grid
    
    # The resolution of the grid
    h = math.sqrt((x[0]-xj[0])**2+(x[1]-xj[1])**2)
    # h = 0.5
    
    # dist = np.empty([1],float)
    dist = math.sqrt((xk[0]-xj[0])**2+(xk[1]-xj[1])**2)
    
    # P = np.empty([2],float)
    # invP = np.empty([2],float) 
    
    if (dist < 0.01*h): #It is the same node
        T = Tj + h*math.sqrt(Q)
        cdir = x - xj
        cdirX = cdir[0]
        cdirY = cdir[1]
    else:
        # Be careful, here Q is cost**2
        dTjk = np.abs(Tj - Tk)
        sqValue = dTjk**2 + 3*h**2*Q - 4*dTjk
        if sqValue < 0.0:
            T = min(Tj,Tk) +  h/math.sqrt(Q)
        else:  
            T = min(Tj,Tk) +  (dTjk + math.sqrt(sqValue))/2
        # T = (Tj + Tk + math.sqrt((Tj+Tk)**2 + 3*h**2*Q - 4*(Tj**2+Tk**2-Tj*Tk)))/2
        
        # This is (x - xj)/h = [djx, djy]
        djx = (x[0] - xj[0])/h
        djy = (x[1] - xj[1])/h
        dkx = (x[0] - xk[0])/h
        dky = (x[1] - xk[1])/h
        
        # The differences in total cost with respect to the neighbors
        dTj = (T - Tj)/h
        dTk = (T - Tk)/h
        
        # Given P = [djk, djy; dkx, dky]
        detP = djx*dky - djy*dkx
        
        # Characteristic Direction = P^-1 dot [dTj, dTk]
        
        # P = np.array([(x - xj)/h, (x - xk)/h])
        # invP = np.linalg.inv(P)
        # Tarray = np.array([dTj,dTk])
        # cdir = invP.dot(Tarray)
        
        cdirX = dky/detP*dTj - djy/detP*dTk
        cdirY = -dkx/detP*dTj + djx/detP*dTk
        
    # dirCharacteristic = math.atan2(cdir[1],cdir[0])
    dirCharacteristic = math.atan2(cdirY,cdirX)
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

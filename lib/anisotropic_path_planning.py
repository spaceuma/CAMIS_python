# -*- coding: utf-8 -*-
"""
Author: J. Ricardo Sanchez Ibanez

Python Implementation of bi-OUM

"""

import matplotlib.pyplot as plt
import numpy as np
import operator
import time
import bisect
import lib.hexgrid as hexgrid
import lib.camis as camis

def computeTmap(VCmap,aspectMap,anisotropyMap,goal,start,Xmap,Ymap,res):
    VCmapG = np.ones_like(VCmap)
    VCmapG[0] = VCmap[0]
    VCmapG[1] = VCmap[1]
    VCmapG[2] = -VCmap[2]
    VCmapG[3] = -VCmap[3]
    
    # State Maps
    #  State -1 = far
    #  State 0 = narrow
    #  State 1 = accepted
    #  State 2 = closed
    stateMapG = -1*np.ones_like(anisotropyMap);
    stateMapS = stateMapG;
    # Define nodeTargets
    nodeTargetG = [goal[0],goal[1]]
    nodeTargetS = [start[0],start[1]]
    
    stateMapG[nodeTargetG[1],nodeTargetG[0]] = 1
    stateMapS[nodeTargetS[1],nodeTargetS[0]] = 1
    
    # Both closedMaps are created
    stateMapG = np.ones_like(anisotropyMap)*(-1)
    stateMapG[nodeTargetG[1],nodeTargetG[0]] = 1
    stateMapG[np.where(anisotropyMap == np.inf)] = 2
    stateMapG[np.where(anisotropyMap == np.nan)] = 2
    stateMapS = np.ones_like(anisotropyMap)*(-1)
    stateMapS[nodeTargetS[1],nodeTargetS[0]] = 1
    stateMapS[np.where(anisotropyMap == np.inf)] = 2
    stateMapS[np.where(anisotropyMap == np.nan)] = 2
    
    # Direction Maps are initialized
    dirMapG = np.ones_like(anisotropyMap)*np.inf
    dirMapS = np.ones_like(anisotropyMap)*np.inf
    
    # Narrow Bands are initialized
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
    NClist, stateMapG = getNewConsidered(nodeTargetG,stateMapG)   
    TmapG, nbTG, nbNodesG= updateNode(NClist, nbTG, nbNodesG, dirMapG, TmapG, stateMapG, VCmapG, aspectMap, anisotropyMap,Xmap,Ymap,res,goal)
    
    NClist, stateMapS = getNewConsidered(nodeTargetS,stateMapS)
    TmapS, nbTS, nbNodesS = updateNode(NClist, nbTS, nbNodesS, dirMapS, TmapS, stateMapS, VCmap, aspectMap, anisotropyMap,Xmap,Ymap,res,start)

    iter = 1
    size = np.size(anisotropyMap)
    
#    fig, axes = plt.subplots(constrained_layout=True)
#    cc = axes.contourf(Xmap, Ymap, stateMapG, 100)
#    fig.colorbar(cc,location='bottom')
#    axes.set_aspect('equal')
#    plt.show()
#    
    
    while nbNodesG or nbNodesS:
        if nbNodesG:
            nodeTargetG, nbTG, nbNodesG = getMinNB(nbTG, nbNodesG)
            stateMapG[nodeTargetG[1],nodeTargetG[0]] = 1
            NClist, stateMapG = getNewConsidered(nodeTargetG,stateMapG)
            stateMapG = updateStateMap(nodeTargetG,stateMapG)
            TmapG, nbTG, nbNodesG = updateNode(NClist, nbTG, nbNodesG, dirMapG, TmapG, stateMapG, VCmapG, aspectMap, anisotropyMap,Xmap,Ymap,res)
        if nbNodesS:
            nodeTargetS, nbTS, nbNodesS = getMinNB(nbTS, nbNodesS)
            stateMapS[nodeTargetS[1],nodeTargetS[0]] = 1
            NClist, stateMapS = getNewConsidered(nodeTargetS,stateMapS)
            stateMapS = updateStateMap(nodeTargetS,stateMapS)
            TmapS, nbTS, nbNodesS = updateNode(NClist, nbTS, nbNodesS, dirMapS, TmapS, stateMapS, VCmap, aspectMap, anisotropyMap,Xmap,Ymap,res)
        if stateMapS[nodeTargetG[1],nodeTargetG[0]] == 2:
            nodeLink = nodeTargetG
            break
        if stateMapG[nodeTargetS[1],nodeTargetS[0]] == 2:
            nodeLink = nodeTargetS
            break
#        iter = iter + 2
#        print('Completed: ' + "{0:.2f}".format(100*iter/size) + ' %')
#        print(iter)
        
    return TmapG, TmapS, dirMapG, dirMapS, nodeLink,stateMapG,stateMapS

#def initialUpdateNode(nodeTarget, costMap, Tmap, nbT, nbNodes, closedMap):
    

def updateStateMap(nodeTarget,stateMap):
    NN = getNeighbours(nodeTarget)
    NN.append(nodeTarget)
    accepted = getAccepted(NN,stateMap)
    N = []
    for i in range(len(accepted)):
        N = getConsidered(getNeighbours(accepted[i]),stateMap)
        if len(N) == 0:
            stateMap[accepted[i][1],accepted[i][0]] = 2
    return stateMap

def getAccepted(nodeList, stateMap):
    acceptedList = []
    for i in range(len(nodeList)):  
        if stateMap[nodeList[i][1],nodeList[i][0]] == 1:
            acceptedList.append(nodeList[i])
    return acceptedList

def getConsidered(nodeList, stateMap):
    consideredList = []
    for i in range(len(nodeList)):     
        if stateMap[nodeList[i][1],nodeList[i][0]] == 0:
            consideredList.append(nodeList[i])
    return consideredList

def getNewConsidered(nodeTarget,stateMap):
    NClist = []
    for i in range(1,7):
        if i == 1:
            nodeChild = np.add(nodeTarget,[0,-1])
        elif i == 2:
            nodeChild = np.add(nodeTarget,[0,1])
        elif i == 3:
            nodeChild = np.add(nodeTarget,[-1,0])
        elif i == 4:
            nodeChild = np.add(nodeTarget,[1,0])
        elif i == 5:
            nodeChild = np.add(nodeTarget,[1,-1])
        elif i == 6:
            nodeChild = np.add(nodeTarget,[-1,1])
        if stateMap[nodeChild[1],nodeChild[0]] == -1:
            stateMap[nodeChild[1],nodeChild[0]] = 0
            NClist.append(nodeChild)
    return NClist, stateMap



def getEikonal(Thor,Tver,cost):
    if Thor == np.inf:
        if Tver == np.inf:
            return np.inf
        else:
            return Tver+cost
    if Tver == np.inf:
        return Thor+cost
    else:
        if cost < np.abs(Thor-Tver):
            return np.minimum(Thor,Tver)+cost
        else:
            return .5*(Thor+Tver+np.sqrt(2*np.power(cost,2) - np.power(Thor-Tver,2)))

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

def updateNode(NClist, nbT, nbNodes, dirMap, Tmap, stateMap, VCmap, aspectMap, anisotropyMap,Xmap,Ymap,res, startingNode = []):
    #nList = getNeighbours(nodeTarget, closedMap)
    for i in range(len(NClist)):
        nodeTarget = NClist[i]
        aspect = aspectMap[:,nodeTarget[1],nodeTarget[0]]
        anisotropy = anisotropyMap[nodeTarget[1],nodeTarget[0]]
        nfPairs = []
        if len(startingNode)!= 0:
            nfPairs = np.concatenate((startingNode,startingNode))
        else:
            R = int(np.ceil(anisotropy) + 1)
            afList = []
            for j in range(-R,R+1):
                for k in range(-R,R+1):
                    try:
                        if stateMap[nodeTarget[1]+j,nodeTarget[0]+k]==1 and \
                        computeDistance(nodeTarget+[k,j],nodeTarget,Xmap,Ymap) <= 2*res*anisotropy:
                            afList.append([nodeTarget[0]+k,nodeTarget[1]+j])
                    except:
                        pass
            localAFPairs = []
            SS = []
            SS[:] = afList
            while (len(SS)!=0):
                ss = SS[0]
                del SS[0]
                for j in range(len(SS)):
                    if (computeDistance(ss,SS[j],Xmap,Ymap) <= 2*res + np.finfo(float).eps):
                        localAFPairs.append(np.concatenate((ss,SS[j])))
            nfPairs = []
            for j in range(len(localAFPairs)):
                if checkNF(localAFPairs[j],nodeTarget,anisotropy,res):
                    nfPairs.append(localAFPairs[j])
            if len(nfPairs) == 0:
                for j in range(len(afList)):
                    nfPairs.append(np.concatenate((afList[j],afList[j])))
        Q1 = VCmap[0,nodeTarget[1],nodeTarget[0]]
        Q2 = VCmap[1,nodeTarget[1],nodeTarget[0]]
        D1 = VCmap[2,nodeTarget[1],nodeTarget[0]]
        D2 = VCmap[3,nodeTarget[1],nodeTarget[0]]
        T,direc = computeT(nodeTarget, nfPairs, Q1, Q2, D1, D2, aspect, Tmap, dirMap,Xmap,Ymap)
        if np.isinf(Tmap[nodeTarget[1],nodeTarget[0]]):
                nIndex = bisect.bisect_left(nbT,T)
                nbT.insert(nIndex,T)
                nbNodes.insert(nIndex, nodeTarget)
                Tmap[nodeTarget[1],nodeTarget[0]] = T
                dirMap[nodeTarget[1],nodeTarget[0]] = direc
        else:
            if T < Tmap[nodeTarget[1],nodeTarget[0]]:
                tempT = Tmap[nodeTarget[1],nodeTarget[0]]
                nIndex = bisect.bisect_left(nbT,tempT)
                nIndex = next(x for x,n in enumerate(nbNodes[nIndex-1:],nIndex) if np.array_equal(nodeTarget,nbNodes[x]))
                del nbT[nIndex]
                del nbNodes[nIndex]
                nIndex = bisect.bisect_left(nbT,T)
                nbT.insert(nIndex,T)
                nbNodes.insert(nIndex, nodeTarget)
                Tmap[nodeTarget[1],nodeTarget[0]] = T  
                dirMap[nodeTarget[1],nodeTarget[0]] = direc
    return Tmap, nbT, nbNodes

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
        return 'false'
    elif p0 <= 0 or p1 <= 0:
        return 'true'
    elif amin >= 0 and amin <= 1:
        return 'true'
    else:
        return 'false'

def computeDistance(nodeA, nodeB, Xmap, Ymap):
    dx = Xmap[nodeA[1],nodeA[0]] - Xmap[nodeB[1],nodeB[0]]
    dy = Ymap[nodeA[1],nodeA[0]] - Ymap[nodeB[1],nodeB[0]]
    return np.sqrt(dx**2+dy**2)

def computeT(nodeTarget, nfPairs, Q1, Q2, D1, D2, aspect, Tmap, dirMap,Xmap,Ymap):
    if np.isnan(aspect[0]) or np.isnan(aspect[1]):
        aspect = [0,0]
    T = np.inf
    direc = np.nan
    x = np.asarray([Xmap[nodeTarget[1],nodeTarget[0]],Ymap[nodeTarget[1],nodeTarget[0]]])
    if len(np.shape(nfPairs)) == 1:
        xj = np.asarray([Xmap[nfPairs[1],nfPairs[0]],Ymap[nfPairs[1],nfPairs[0]]])
        xk = np.asarray([Xmap[nfPairs[3],nfPairs[2]],Ymap[nfPairs[3],nfPairs[2]]])
        Tj = Tmap[nfPairs[1],nfPairs[0]]
        Tk = Tmap[nfPairs[3],nfPairs[2]]
        preT, preDir = optimizeCost(x,xj,xk,Tj,Tk,Q1,Q2,D1,D2,aspect)
        if preT < T:
            T = preT
            direc = preDir
    else:
        for i in range(len(nfPairs)):
            xj = np.asarray([Xmap[nfPairs[i][1],nfPairs[i][0]],Ymap[nfPairs[i][1],nfPairs[i][0]]])
            xk = np.asarray([Xmap[nfPairs[i][3],nfPairs[i][2]],Ymap[nfPairs[i][3],nfPairs[i][2]]])
            Tj = Tmap[nfPairs[i][1],nfPairs[i][0]]
            Tk = Tmap[nfPairs[i][3],nfPairs[i][2]]
            preT, preDir = optimizeCost(x,xj,xk,Tj,Tk,Q1,Q2,D1,D2,aspect)
            if preT < T:
                T = preT
                direc = preDir
    return T,direc

def optimizeCost(x,xj,xk,Tj,Tk,Q1,Q2,D1,D2,aspect):
    a = 0
    b = 1
    nIter = 30
    tau=(np.sqrt(5)-1)/2;
    epsilon1=a+(1-tau)*(b-a)
    epsilon2=a+tau*(b-a)
    beta1 = camis.computeBeta(aspect,x-epsilon1*xj-(1-epsilon1)*xk)
    f_x1= camis.getCAMIScost(beta1,Q1,Q2,D1,D2) + epsilon1*Tj + (1-epsilon1)*Tk
    beta2 = camis.computeBeta(aspect,x-epsilon2*xj-(1-epsilon2)*xk)
    f_x2= camis.getCAMIScost(beta2,Q1,Q2,D1,D2) + epsilon2*Tj + (1-epsilon2)*Tk
    accuracy = 0.1
    k = 0
    while np.abs(b-a)>accuracy and k<nIter:
        k = k+1
        if f_x1 < f_x2:
            b = epsilon2
            epsilon2 = epsilon1
            epsilon1 = a+(1-tau)*(b-a)
        else:
            a = epsilon1
            epsilon1 = epsilon2
            epsilon2 = a+tau*(b-a)
        beta1 = camis.computeBeta(aspect,x-epsilon1*xj-(1-epsilon1)*xk)
        f_x1= camis.getCAMIScost(beta1,Q1,Q2,D1,D2) + epsilon1*Tj + (1-epsilon1)*Tk
        beta2 = camis.computeBeta(aspect,x-epsilon2*xj-(1-epsilon2)*xk)
        f_x2= camis.getCAMIScost(beta2,Q1,Q2,D1,D2) + epsilon2*Tj + (1-epsilon2)*Tk
    if f_x1 < f_x2:
        minEpsilon = epsilon1
        minT = f_x1
    else:
        minEpsilon = epsilon2
        minT = f_x2
    
    minDirVector = (x-minEpsilon*xj-(1-minEpsilon)*xk)
    minDirVector = minDirVector/np.linalg.norm(minDirVector)
    minDir = np.arctan2(minDirVector[1],minDirVector[0])
    return minT,minDir

    
def getMinNB(nbT,nbNodes):
    nodeTarget = nbNodes.pop(0)
    del nbT[0]
    return nodeTarget, nbT, nbNodes

# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 16:56:51 2019

@author: rsanchez
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import lib.camislib as camislib
import scipy.interpolate as interp
import lib.anisotropic_path_planning as ap

from scipy import signal
from scipy import ndimage
from time import time
from matplotlib import cbook
from matplotlib import cm
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

rad2deg = 180/np.pi


# Proccessed DEM
class PDEM:
    def __init__(self, DEM, demRes, planRes, offset):
        self.elevationMap = DEM
        self.size = DEM.shape
        self.xMap, self.yMap = \
          np.meshgrid(np.linspace(offset[0],self.size[0]-1,self.size[0]), \
                      np.linspace(offset[1],self.size[1]-1,self.size[1]))
        self.demRes = demRes
        self.planRes = planRes
        self.offset = offset
        self.computeSlope()
        self.computeLaplacian()
        self.maxSlope = np.max(self.slopeMap)
        self.xMin = self.xMap[0][0]
        self.yMin = self.yMap[0][0]
    def smoothMap(self, radius):
        r = int(radius/self.demRes)
        r = r + 1 - r%2
        y,x = np.ogrid[-r: r+1, -r: r+1]
        convMatrix = x**2+y**2 <= r**2
        convMatrix = convMatrix.astype(float)
        convMatrix = convMatrix/convMatrix.sum()
        smoothDEM = signal.medfilt2d(self.elevationMap,r)
        smoothDEM = signal.convolve2d(self.elevationMap, convMatrix, \
                                      mode='same', boundary='symm')
        self.elevationMap = signal.convolve2d(smoothDEM, convMatrix, \
                                      mode='same', boundary='symm')
        
        self.computeSlope()
        
    def computeSlope(self):
        self.gradientY, self.gradientX = np.gradient(self.elevationMap,\
                                                     self.demRes,\
                                                     self.demRes)
        self.gradientMap = np.sqrt(self.gradientX**2+self.gradientY**2)
        self.slopeMap = np.arctan(self.gradientMap)
        self.aspectMap = np.arctan2(self.gradientY,self.gradientX)+np.pi
        self.aspectX = np.cos(self.aspectMap)
        self.aspectY = np.sin(self.aspectMap)
        
    def computeLaplacian(self):
        self.laplacianY, self.laplacianX = np.gradient(self.gradientMap,\
                                                     self.demRes,\
                                                     self.demRes)
        self.laplacianMap = self.laplacianX + self.laplacianY
        
    def show3dDEM(self):
        fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
        ls = LightSource(270, 45)
        # To use a custom hillshading mode, override the built-in shading and pass
        # in the rgb colors of the shaded surface calculated from "shade".
        rgb = ls.shade(self.elevationMap, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
        ax.plot_surface(self.xMap, self.yMap, self.elevationMap, rstride=1, cstride=1, facecolors=rgb,
                               linewidth=0, antialiased=False, shade=False)
        ax.auto_scale_xyz([0, 51], [0, 51], [0, 51])
        plt.show()
    def showMap(self, opt, fig, axes):
        if   opt == 'elevation':
            cc = axes.contourf(self.elevationMap, 100, cmap = cm.gist_earth)
        elif opt == 'slope-rad':
            cc = axes.contourf(self.slopeMap, 100, cmap = 'plasma')
        elif opt == 'slope-deg':
            cc = axes.contourf(rad2deg*self.slopeMap, 100, cmap = 'plasma')
        elif opt == 'aspect-rad':
            cc = axes.contourf(self.aspectMap, 20,cmap = 'hsv')
            cc.set_clim(0,2*np.pi)
        elif opt == 'aspect-deg':
            cc = axes.contourf(rad2deg*self.aspectMap, 20,cmap = 'hsv')
            cc.set_clim(0,360)
        elif opt == 'laplacian':
            cc = axes.contourf(self.laplacianMap, 100, cmap = 'plasma')
        elif opt == 'laplacian-abs':
            cc = axes.contourf(np.abs(self.laplacianMap), 100, cmap = 'plasma')
        cbar = fig.colorbar(cc,location='bottom')
        cbar.set_label('Height (m)')
        axes.set_aspect('equal')
        plt.show()



class AnisotropicMap(PDEM):
    @classmethod
    def from_pdem(cls, parentPdem, camis):
        return cls(parentPdem.elevationMap, parentPdem.demRes,\
                   parentPdem.planRes, parentPdem.offset, camis)
    
    def __init__(self, DEM, demRes, planRes, offset):
        super().__init__(DEM, demRes, planRes, offset)
        self.computeHexGrid()
    
    def computeHexGrid(self):
        DX = self.xMap[-1,-1] - self.xMap[0,0]
        DY = self.yMap[-1,-1] - self.yMap[0,0]
        JMax = math.ceil(2*DY/(math.sqrt(3)*self.planRes))
        IMax = math.ceil((DX+DY/math.sqrt(3))/self.planRes)
        II,JJ = np.meshgrid(np.linspace(0,IMax,IMax+1),np.linspace(0,JMax,JMax+1))
        self.hexXmap = self.xMap[0,0] + self.planRes*(II + .5*JJ) - DY/math.sqrt(3)
        self.hexYmap = self.yMap[0,0] + self.planRes*math.sqrt(3)/2*JJ
        XX,YY = np.meshgrid(np.linspace(0,np.ceil(DX),np.ceil(DX)+1),np.linspace(0,np.ceil(DY),np.ceil(DY)+1))
        self.xy2J = 2*YY/(np.sqrt(3)*self.planRes)
        self.xy2I = (DY/np.sqrt(3)+ XX)/self.planRes-0.5*self.xy2J
    
    def computeVecCostMap(self, camis):
        self.camis = camis
        obstacleMap = np.zeros_like(self.slopeMap)
        obstacleMap[0,:] = 1
        obstacleMap[-1,:] = 1
        obstacleMap[:,0] = 1
        obstacleMap[:,-1] = 1
        proximityMap = ndimage.morphology.distance_transform_edt(1-obstacleMap)
        proximityMap[np.where(proximityMap[:]<0)] = 0
        self.obstacleMap = obstacleMap
        self.proximityMap = proximityMap
        
        points = np.zeros((self.xMap.size,2))
        points[:,0] = self.xMap.flatten()
        points[:,1] = self.yMap.flatten()
        hexSlopeMap = interp.griddata(points, self.slopeMap.flatten(),\
                                      (self.hexXmap, self.hexYmap), \
                                      method='nearest')
        hexSlopeMap[np.where(self.hexXmap < self.xMap[0,0])] = np.nan
        hexSlopeMap[np.where(self.hexXmap > self.xMap[-1,-1])] = np.nan
        hexSlopeMap[np.where(self.hexYmap < self.yMap[0,0])] = np.nan
        hexSlopeMap[np.where(self.hexYmap > self.yMap[-1,-1])] = np.nan
        
        init = time()


        cdRoots = self.camis.cdRoots
        caRoots = self.camis.caRoots
        cl1Roots = self.camis.cl1Roots
        cl2Roots = self.camis.cl2Roots
        anicoLUT = self.camis.anicoLUT
        vectorialData = camislib.getVectorialCostMap(rad2deg*hexSlopeMap,cdRoots,caRoots,cl1Roots,\
                                          cl2Roots,anicoLUT)
        
        
        AspectMap = np.zeros([2,hexSlopeMap.shape[0],hexSlopeMap.shape[1]])
        AspectMap[0] = interp.griddata(points, self.aspectX.flatten(),\
                 (self.hexXmap, self.hexYmap), method='nearest')
        AspectMap[1] = interp.griddata(points, self.aspectY.flatten(),\
                 (self.hexXmap, self.hexYmap), method='nearest')
        ProximityMap = interp.griddata(points, self.proximityMap.flatten(),\
                                       (self.hexXmap, self.hexYmap), method='nearest')
        AspectMap[0][np.where(np.isnan(hexSlopeMap))] = np.nan
        AspectMap[1][np.where(np.isnan(hexSlopeMap))] = np.nan
        ProximityMap[np.where(np.isnan(hexSlopeMap))] = np.nan
        obstacleMask = ProximityMap == 0
        
        AnisotropyMap = vectorialData[0][:][:]
        AnisotropyMap[obstacleMask] = np.inf
        
        VCMap = np.zeros([4,AnisotropyMap.shape[0],AnisotropyMap.shape[1]])
        
        Q1 = vectorialData[1][:][:]
        Q1[obstacleMask] = np.inf
        Q2 = vectorialData[2][:][:]
        Q2[obstacleMask] = np.inf
        D1 = vectorialData[3][:][:]
        D1[obstacleMask] = np.inf
        D2 = vectorialData[4][:][:]
        D2[obstacleMask] = np.inf
        
        VCMap[0] = Q1
        VCMap[1] = Q2
        VCMap[2] = D1
        VCMap[3] = D2
        
        self.VCMap = VCMap
        self.hexSlopeMap = hexSlopeMap
        self.hexAspectMap = AspectMap
        self.hexAnisotropyMap = AnisotropyMap
        
        print('Elapsed time to compute the Vectorial Cost Map: '+str(time()-init))
    
    def executePlanning(self, goal, start):
        init = time()
        IJ2XY = np.zeros([2,self.hexXmap.shape[0],self.hexYmap.shape[1]])
        IJ2XY[0] = self.hexXmap
        IJ2XY[1] = self.hexYmap
        XY2IJ = np.zeros([2,self.xy2I.shape[0],self.xy2J.shape[1]])
        XY2IJ[0] = self.xy2I
        XY2IJ[1] = self.xy2J
        ijStart = np.round(XY2IJ[:,start[1],start[0]]).astype(int)
        ijGoal = np.round(XY2IJ[:,goal[1],goal[0]]).astype(int)
        self.Tmap, dirMap, stateMap, maxAnisoMap = \
        ap.computeTmap(self.VCMap, self.hexAspectMap, self.hexAnisotropyMap,\
                       ijGoal, ijStart, self.hexXmap, self.hexYmap, self.planRes)
        print('Elapsed time to compute the Total Cost Map: '+str(time()-init))
        
        
#        startWaypoint = IJ2XY[:,start[1],start[0]]
#        goalWaypoint = IJ2XY[:,goal[1],goal[0]]
        self.IJ2XY = IJ2XY
        self.XY2IJ = XY2IJ
        path,uu = ap.getPath(dirMap, IJ2XY, XY2IJ, start, goal, self.xMin, self.yMin, self.planRes)
        self.path = np.asarray(path)
        self.getPathData()
        self.dirMap = dirMap
        
    def getPathData(self):
        pathElevation = []
        pathSlope = []
        pathAspectX = []
        pathAspectY = []
        pathTravDist = []
        pathPitch = []
        pathRoll = []
        pathHeading = []
        pathAspect = []
        pathCost = []
        pathSegment = []
        pathEstimatedTotalCost = []
        pathComputedTotalCost = []
        for index, waypoint in enumerate(self.path):
            pathElevation.append(ap.interpolatePoint(waypoint,self.elevationMap))
#            pathSlope.append(ap.interpolatePoint(waypoint,self.slopeMap))
#            pathAspect.append(np.arctan2(ap.interpolatePoint(waypoint,self.aspectY), ap.interpolatePoint(waypoint,self.aspectX)))
            pathSlope.append(ap.getTriInterpolation(waypoint,self.hexSlopeMap,self.XY2IJ,self.IJ2XY,self.planRes, self.xMin, self.yMin))
            pathAspectX.append(ap.getTriInterpolation(waypoint,self.hexAspectMap[0],self.XY2IJ,self.IJ2XY,self.planRes, self.xMin, self.yMin))
            pathAspectY.append(ap.getTriInterpolation(waypoint,self.hexAspectMap[1],self.XY2IJ,self.IJ2XY,self.planRes, self.xMin, self.yMin))
            pathAspect.append(np.arctan2(pathAspectY[index],pathAspectX[index]))
            if index == 0:
                pathTravDist.append(0)
            else:
                A = [self.path[index][0], self.path[index][1], pathElevation[index]]
                A = np.asarray(A)
                B = [self.path[index-1][0], self.path[index-1][1], pathElevation[index-1]]
                B = np.asarray(B)
                pathTravDist.append(pathTravDist[index-1] + np.linalg.norm(A-B))
        for index, waypoint in enumerate(self.path):
            if index == self.path.shape[0]-1:
                pathPitch.append(np.nan) #Fix this!
                pathHeading.append(pathHeading[-1])
                pathSegment.append(0)
            elif index == 0:
                pathPitch.append(np.nan)
                A = self.path[index+1] - self.path[index]
                pathSegment.append(np.linalg.norm(A))
                pathHeading.append(np.arctan2(A[1],A[0]))
            else:
                pathPitch.append((pathElevation[index]-pathElevation[index+1])/(pathTravDist[index+1]-pathTravDist[index]))
                A = self.path[index+1] - self.path[index]
                pathSegment.append(np.linalg.norm(A))
                pathHeading.append(np.arctan2(A[1],A[0]))
            pathCost.append(self.camis.getCost(rad2deg*pathSlope[index],pathAspect[index],pathHeading[index]))
#        N = 5
#        pathPitch = np.convolve(pathPitch, np.ones((N,))/N, mode='same')
        for index, waypoint in enumerate(self.path):
            pathRoll.append(np.arccos(np.cos(pathSlope[index])/np.cos(pathPitch[index])))
            if np.isnan(pathRoll[index]): #This can happen due to numerical errors
                pathRoll[index] = 0.0
            pathEstimatedTotalCost.append(ap.getTriInterpolation(waypoint,self.Tmap,self.XY2IJ,self.IJ2XY,self.planRes, self.xMin, self.yMin))
            if index == 0:
                pathComputedTotalCost.append(0)
            else:
                pathComputedTotalCost.append(pathComputedTotalCost[index-1]+0.5*(pathCost[index-1]+pathCost[index])*pathSegment[index-1])
            
        self.pathElevation = pathElevation
        self.pathSlope = pathSlope
        self.pathTravDist = pathTravDist
        self.pathPitch = pathPitch
        self.pathRoll = pathRoll
        self.pathAspect = pathAspect
        self.pathHeading = pathHeading
        self.pathCost = pathCost
        self.pathSegment = pathSegment
        self.pathComputedTotalCost = pathComputedTotalCost
        self.pathEstimatedTotalCost = pathEstimatedTotalCost
        
    def showVecCostMap(self, index):
        fig, ax = plt.subplots()
        cc = ax.contourf(self.hexXmap, self.hexYmap, self.VCMap[index], 100, cmap = 'plasma')
        fig.colorbar(cc)
        ax.set_aspect('equal')
        plt.show()
        
    def showHexSlopeMap(self):
        fig, ax = plt.subplots()
        cc = ax.contourf(self.hexXmap, self.hexYmap, rad2deg*self.hexSlopeMap, 100, cmap = 'plasma')
        fig.colorbar(cc)
        ax.set_aspect('equal')
        plt.show()
        
    def showResults(self):
        fig, ax = plt.subplots()
        cc = ax.contourf(self.hexXmap, self.hexYmap, self.Tmap, 100, cmap = 'nipy_spectral', alpha = .5)
        ax.contour(self.hexXmap, self.hexYmap, self.Tmap, 100, cmap = 'nipy_spectral')
        ax.quiver(self.hexXmap, self.hexYmap,np.cos(self.dirMap), np.sin(self.dirMap))
        ax.plot(self.path[:,0],self.path[:,1],'r')
        cbar = fig.colorbar(cc)
        cbar.set_label('Total Cost')
        ax.set_aspect('equal')
        plt.show()
        
        fig, ax = plt.subplots()
        ax.plot(self.pathTravDist, [x*rad2deg for x in self.pathSlope], linestyle='solid')
        ax.plot(self.pathTravDist, [x*rad2deg for x in self.pathPitch], linestyle='dotted')
        ax.plot(self.pathTravDist, [x*rad2deg for x in self.pathRoll], linestyle='dashed')
        ax.set_aspect('equal')
        plt.show()
        
        fig, ax = plt.subplots()
        ax.plot(self.pathTravDist, self.pathComputedTotalCost, linestyle='solid')
        ax.plot(self.pathTravDist, [self.pathEstimatedTotalCost[0]-x for x in self.pathEstimatedTotalCost], linestyle='dotted')
#        ax.set_aspect('equal')
        ax.legend(('Computed', 'Estimated'))
        plt.show()
        
        
    def showPath(self, fig, axes, color):
        axes.plot(self.path[:,0],self.path[:,1], color)
        plt.show()
        
    def showPathData(self, opt, fig, axes, color):
        if   opt == 'elevation':
            axes.plot(self.pathTravDist, self.pathElevation, color)
        elif opt == 'slope':
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathSlope], color)
        elif opt == 'pitch':
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathPitch], color)
        elif opt == 'roll':
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathRoll], color)
        elif opt == 'cost':
            axes.plot(self.pathTravDist, self.pathCost, color)
        elif opt == 'total-cost':
            axes.plot(self.pathTravDist, self.pathComputedTotalCost, color)
        elif opt == 'full-orientation':
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathSlope], color, linestyle='solid')
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathPitch], color, linestyle='dotted')
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathRoll], color, linestyle='dashed')
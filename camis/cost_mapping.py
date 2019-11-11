# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 16:56:51 2019

@author: rsanchez
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from .cost_function import getVectorialCostMap
import scipy.interpolate as interp
import camis.anisotropic_path_planning as ap
import matplotlib

from mayavi import mlab

from scipy import signal
from scipy import ndimage
from time import time
from matplotlib import cbook
from matplotlib import cm
from matplotlib.colors import LightSource
from mpl_toolkits.mplot3d import Axes3D

rad2deg = 180/np.pi


# Proccessed DEM
class PDEM:
    def __init__(self, DEM, demRes, planRes, offset):
        self.elevationMap = DEM
        self.size = DEM.shape
        xMap, yMap = \
          np.meshgrid(np.linspace(0,self.size[1]-1,self.size[1]), \
                      np.linspace(0,self.size[0]-1,self.size[0]))
        self.xMap = xMap*demRes
        self.yMap = yMap*demRes
        self.demRes = demRes
        self.planRes = planRes
        self.offset = offset
        self.computeSlope()
        self.computeLaplacian()
        self.maxSlope = np.max(self.slopeMap)
        self.xMin = self.xMap[0][0]
        self.yMin = self.yMap[0][0]
        
    def smoothMap(self, radius):
        self.oldElevationMap = self.elevationMap
        self.oldSlopeMap = self.slopeMap
        self.computeRoughness(radius)
        r = int(radius/self.demRes)
        r = r + 1 - r%2
        y,x = np.ogrid[-r: r+1, -r: r+1]
        convMatrix = x**2+y**2 <= r**2
        convMatrix = convMatrix.astype(float)
#        convMatrix = convMatrix/convMatrix.sum()
#        smoothDEM = signal.medfilt2d(self.elevationMap,r)
        medDEM = ndimage.median_filter(self.elevationMap, footprint=convMatrix, mode='nearest')
#        smoothDEM = signal.convolve2d(smoothDEM, convMatrix, \
#                                      mode='same', boundary='symm')
#        self.computeSlope()
        self.maxSlopeMap = ndimage.maximum_filter(self.slopeMap, footprint=convMatrix, mode='nearest')
        
#        smoothDEM = signal.convolve2d(self.elevationMap, convMatrix/convMatrix.sum(), \
#                                      mode='same', boundary='symm')
#        self.elevationMap = medDEM + (medDEM - smoothDEM)
        
        
        self.elevationMap = signal.convolve2d(self.elevationMap, convMatrix/convMatrix.sum(), \
                                      mode='same', boundary='symm')
        self.computeSlope()
    
    def processSlope(self, sdThreshold):
#        self.computeLaplacian()
        AA = np.ones_like(self.maxSlopeMap)*30.0*np.pi/180.0
        self.slopeMap[np.where(self.sdMap > sdThreshold)] = \
                        np.maximum(AA[np.where(self.sdMap > sdThreshold)],\
                        self.maxSlopeMap[np.where(self.sdMap > sdThreshold)])
    
    def computeSlope(self):
        self.gradientY, self.gradientX = np.gradient(self.elevationMap,\
                                                     self.demRes,\
                                                     self.demRes)
        self.gradientMap = np.sqrt(self.gradientX**2+self.gradientY**2)
        self.slopeMap = np.arctan(self.gradientMap)
        self.aspectMap = np.arctan2(self.gradientY,self.gradientX)+np.pi
        self.aspectX = np.cos(self.aspectMap)
        self.aspectY = np.sin(self.aspectMap)
    
    def computeRoughness(self,radius):
#        self.normalX = self.aspectX/(1+(1/np.tan(self.slopeMap))**2)
#        self.normalY = self.aspectY/(1+(1/np.tan(self.slopeMap))**2)
#        self.normalZ = 1/(1+np.tan(self.slopeMap)**2)
        
        self.normalX = np.cos(self.aspectMap)*np.sin(self.slopeMap)
        self.normalY = np.sin(self.aspectMap)*np.sin(self.slopeMap)
        self.normalZ = np.cos(self.slopeMap)
        
        r = int(radius/self.demRes)
        r = r + 1 - r%2
        y,x = np.ogrid[-r: r+1, -r: r+1]
        convMatrix = x**2+y**2 <= r**2
        convMatrix = convMatrix.astype(float)
        
        
        sumNormalX = signal.convolve2d(self.normalX, convMatrix, \
                                      mode='same', boundary='symm')
        sumNormalY = signal.convolve2d(self.normalY, convMatrix, \
                                      mode='same', boundary='symm')
        sumNormalZ = signal.convolve2d(self.normalZ, convMatrix, \
                                      mode='same', boundary='symm')
        
        self.roughnessMap = np.sqrt( sumNormalX**2 + sumNormalY**2 + sumNormalZ**2 )
        #Standard Deviation
        self.sdMap = np.sqrt( -2*np.log(self.roughnessMap/convMatrix.sum()))*180/np.pi
        
        #Dispersion
        self.dispersionMap = (convMatrix.sum()-self.roughnessMap)/(convMatrix.sum()-1)
        
    def computeLaplacian(self):
        self.laplacianY, self.laplacianX = np.gradient(self.gradientMap,\
                                                     self.demRes,\
                                                     self.demRes)
        self.laplacianMap = self.laplacianX + self.laplacianY
    def getBeta(self, xPos, yPos, heading):
        slope = []
        beta = []
        for index, waypoint in enumerate(xPos):
            try:
                aspectX = ap.interpolatePoint([xPos[index]/self.demRes,yPos[index]/self.demRes],self.aspectX)
                aspectY = ap.interpolatePoint([xPos[index]/self.demRes,yPos[index]/self.demRes],self.aspectY)
            except:
                print('ERROR')
            aspect = np.arctan2(aspectY,aspectX)
            slope.append(ap.interpolatePoint([xPos[index]/self.demRes,yPos[index]/self.demRes],self.slopeMap))
            b = np.arccos(np.cos(heading[index])*np.cos(aspect)+np.sin(heading[index])*np.sin(aspect))
            crossDirection = np.sin(heading[index])*np.cos(aspect)-np.cos(heading[index])*np.sin(aspect)
            if crossDirection >= 0.:
                beta.append(b)
            else:
                beta.append(-b)
        return slope,beta
        
    def show3dDEM(self):
#        fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), constrained_layout=True)
#        ls = LightSource(270, 45)
#        # To use a custom hillshading mode, override the built-in shading and pass
#        # in the rgb colors of the shaded surface calculated from "shade".
#        slopeMap = rad2deg*self.slopeMap
#        vmin=slopeMap.min().min()
#        vmax=slopeMap.max().max()
#        zRange = max(-self.xMap[0][0]+self.xMap[-1][-1],-self.yMap[0][0]+self.yMap[-1][-1])
#        rgb = ls.shade(self.elevationMap, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
#        ax.plot_surface(self.xMap, self.yMap, self.elevationMap, rstride=1, cstride=1, facecolors=rgb,
#                               linewidth=0, antialiased=False, shade=False)
#        cset = ax.contourf(self.xMap, self.yMap, slopeMap, zdir='z', offset=self.elevationMap.max().max()-zRange, cmap=cm.gist_rainbow)
#        ax.view_init(elev=35., azim=-120.)
#        
##        norm = matplotlib.colors.Normalize(vmin, vmax)
##        m = cm.ScalarMappable(cmap=plt.cm.gist_earth, norm=norm)
##        m.set_array([])
##        plt.colorbar(m)
#        cbar = fig.colorbar(cset,orientation='horizontal',shrink=0.5)
#        cbar.set_label('Slope (deg)')
##        ax.autoscale()
##        ax.auto_scale_xyz([self.xMap[0][0], self.xMap[-1][-1]], \
##                          [self.yMap[0][0], self.yMap[-1][-1]], \
##                          [self.elevationMap.min().min(), self.elevationMap.max().max()])
#        ax.set_aspect('equal')
#        ax.set_xlim(self.xMap[0][0], self.xMap[-1][-1])
#        ax.set_ylim(self.yMap[0][0], self.yMap[-1][-1])
#        ax.set_zlim(self.elevationMap.max().max()-zRange, self.elevationMap.max().max())
#        
#        plt.show()
        mlab.figure(size=(640, 800))
        mlab.surf(self.oldElevationMap, colormap='gist_earth', warp_scale='auto')
        mlab.show()
    def showMap(self, opt, fig, axes):
        if   opt == 'elevation':
            cc = axes.contourf(self.xMap, self.yMap, self.elevationMap, 100, cmap = cm.gist_earth)
        if   opt == 'old-elevation':
#            levels = np.linspace(47.0,55.0,25)
#            cc = axes.contourf(self.xMap, self.yMap, self.oldElevationMap, levels=levels, cmap = cm.gist_earth, extend='both')
            cc = axes.contourf(self.xMap, self.yMap, self.oldElevationMap, 50, cmap = cm.gist_earth, extend='both')
            axes.contour(self.xMap, self.yMap, self.oldElevationMap, 50, colors = 'k', alpha=.3)
#            cc.set_clim(47.0,55.0)
            cbar = fig.colorbar(cc)
            cbar.set_label('Height (m)')
        if   opt == 'elevation-contour':
            cc = axes.contourf(self.xMap, self.yMap, self.elevationMap, 100, cmap = cm.gist_earth)
            axes.contour(self.xMap, self.yMap, self.elevationMap, 100, colors = 'k', alpha=.3)
        elif opt == 'slope-rad':
            cc = axes.contourf(self.xMap, self.yMap, self.slopeMap, 100, cmap = 'plasma')
        elif opt == 'slope-deg':
            cc = axes.contourf(self.xMap, self.yMap, rad2deg*self.slopeMap, 100, cmap = 'nipy_spectral')
#            axes.quiver(self.xMap, self.yMap,self.aspectX, self.aspectY,scale = 40)
#            cc.set_clim(0,45.0)
        elif opt == 'roughness':
            cc = axes.contourf(self.xMap, self.yMap, self.roughnessMap, 100, cmap = 'nipy_spectral')
        elif opt == 'standard-deviation':
            cc = axes.contourf(self.xMap, self.yMap, self.sdMap, 100, cmap = 'nipy_spectral')
        elif opt == 'dispersion':
            cc = axes.contourf(self.xMap, self.yMap, self.dispersionMap, 100, cmap = 'nipy_spectral')
        elif opt == 'old-slope-deg':
            cc = axes.contourf(self.xMap, self.yMap, rad2deg*self.oldSlopeMap, 100, cmap = 'nipy_spectral')
#            axes.quiver(self.xMap, self.yMap,self.normalX, self.normalY,scale = 40)
        elif opt == 'aspect-rad':
            cc = axes.contourf(self.xMap, self.yMap, self.aspectMap, 20,cmap = 'hsv')
            cc.set_clim(0,2*np.pi)
        elif opt == 'aspect-deg':
            cc = axes.contourf(self.xMap, self.yMap, rad2deg*self.aspectMap, 20,cmap = 'hsv')
            cc.set_clim(0,360)
        elif opt == 'laplacian':
            cc = axes.contourf(self.xMap, self.yMap, self.laplacianMap, 100, cmap = 'plasma')
        elif opt == 'laplacian-abs':
            cc = axes.contourf(self.xMap, self.yMap, np.abs(self.laplacianMap), 100, cmap = 'nipy_spectral')
            cc.set_clim(0,5.0)
        axes.set_aspect('equal')
        plt.show()



class AnisotropicMap(PDEM):
    @classmethod
    def from_pdem(cls, parentPdem, camis):
        return cls(parentPdem.elevationMap, parentPdem.demRes,\
                   parentPdem.planRes, parentPdem.offset, camis)
    
    def __init__(self, DEM, demRes, planRes, offset):
        init = time()
        super().__init__(DEM, demRes, planRes, offset)
        self.computeHexGrid()
        print('Hexagonal Grid created in ' + str(time()-init))
    
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
        
        # We define the forbidden areas
        obstacleMap = np.zeros_like(self.slopeMap)
        obstacleMap[0,:] = 1
        obstacleMap[-1,:] = 1
        obstacleMap[:,0] = 1
        obstacleMap[:,-1] = 1
#        self.slopeMap[np.where(self.maxSlopeMap > np.pi/180.0*camis.slopeThreshold)] = self.maxSlopeMap[np.where(self.maxSlopeMap > np.pi/180.0*camis.slopeThreshold)]
#        obstacleMap[np.where(self.slopeMap > np.pi/180.0*1.5*camis.slopeThreshold)] = 1
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
        vectorialData = getVectorialCostMap(rad2deg*hexSlopeMap,cdRoots,caRoots,cl1Roots,\
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
        D1[obstacleMask] = 0
        D2 = vectorialData[4][:][:]
        D2[obstacleMask] = 0
        
        VCMap[0] = Q1
        VCMap[1] = Q2
        VCMap[2] = D1
        VCMap[3] = D2
        
        self.VCMap = VCMap
        self.hexSlopeMap = hexSlopeMap
        self.hexAspectMap = AspectMap
        self.hexAnisotropyMap = AnisotropyMap
        
        print('Elapsed time to compute the Vectorial Cost Map: '+str(time()-init))
    
    # Executing OUM to compute the path
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
        self.Tmap, dirMap, stateMap = \
        ap.computeTmap(self.VCMap, self.hexAspectMap, self.hexAnisotropyMap,\
                       ijGoal, ijStart, self.hexXmap, self.hexYmap, self.planRes)
        print('Elapsed time to compute the Total Cost Map: '+str(time()-init))
        
        
#        startWaypoint = IJ2XY[:,start[1],start[0]]
#        goalWaypoint = IJ2XY[:,goal[1],goal[0]]
        self.IJ2XY = IJ2XY
        self.XY2IJ = XY2IJ
        self.dirMap = dirMap
        path,uu = ap.getPath(dirMap, IJ2XY, XY2IJ, start, goal, self.xMin, self.yMin, self.planRes)
        self.path = np.asarray(path)
        self.getPathData()
    
    # Executing BiOUM to compute the path faster
    def executeBiPlanning(self, goal, start):
        init = time()
        IJ2XY = np.zeros([2,self.hexXmap.shape[0],self.hexYmap.shape[1]])
        IJ2XY[0] = self.hexXmap
        IJ2XY[1] = self.hexYmap
        XY2IJ = np.zeros([2,self.xy2I.shape[0],self.xy2J.shape[1]])
        XY2IJ[0] = self.xy2I
        XY2IJ[1] = self.xy2J
        ijStart = np.round(XY2IJ[:,start[1],start[0]]).astype(int)
        ijGoal = np.round(XY2IJ[:,goal[1],goal[0]]).astype(int)
        # TmapG, TmapS, dirMapG, dirMapS, nodeLink, stateMapG, stateMapS
        self.TmapG, self.TmapS, self.dirMapG, self.dirMapS, nodeLink, stateMapG, stateMapS, self.dirLinkG, self.dirLinkS = \
        ap.computeBiTmap(self.VCMap, self.hexAspectMap, self.hexAnisotropyMap,\
                       ijGoal, ijStart, self.hexXmap, self.hexYmap, self.planRes)
        print('Elapsed time to compute the Total Cost Map: '+str(time()-init))
        self.linkNode = nodeLink
        linkCoord = IJ2XY[:,nodeLink[1],nodeLink[0]]
        
        self.pathG,uu = ap.getPath(self.dirMapG, IJ2XY, XY2IJ, linkCoord, goal, self.xMin, self.yMin, self.planRes)
        self.pathS,uu = ap.getPath(self.dirMapS, IJ2XY, XY2IJ, linkCoord, start, self.xMin, self.yMin, self.planRes)
        
        self.IJ2XY = IJ2XY
        self.XY2IJ = XY2IJ
        
        TlinkG = self.TmapG[nodeLink[1],nodeLink[0]]
        TlinkS = self.TmapS[nodeLink[1],nodeLink[0]]
        
        self.path = np.concatenate((np.flipud(np.asarray(self.pathS)),np.asarray(self.pathG)))
        self.Tmap = np.zeros_like(self.TmapG)
        self.Tmap[:] = self.TmapG
        self.Tmap[np.where(self.TmapS != np.inf)] = TlinkS + TlinkG - self.TmapS[np.where(self.TmapS != np.inf)]
        self.Tmap[np.where(self.Tmap < 0.0)] = np.inf
        self.getPathData()
        self.linkPos = linkCoord
        
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
        pathBeta = []
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
            pathBeta.append(np.arccos(np.cos(pathHeading[index])*np.cos(pathAspect[index])+np.sin(pathHeading[index])*np.sin(pathAspect[index])))
            if np.isnan(pathRoll[index]): #This can happen due to numerical errors
                pathRoll[index] = 0.0
            if (index == self.path.shape[0] - 1):
                print('stop')
            pathEstimatedTotalCost.append(ap.getTriLowest(waypoint,self.Tmap,self.XY2IJ,self.IJ2XY,self.planRes, self.xMin, self.yMin))
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
        self.pathBeta = pathBeta
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
        cc = ax.contourf(self.hexXmap, self.hexYmap, rad2deg*self.hexSlopeMap, 100, cmap = 'nipy_spectral')
        ax.quiver(self.hexXmap, self.hexYmap,self.hexAspectMap[0], self.hexAspectMap[1],scale = 30)
        fig.colorbar(cc)
        ax.set_aspect('equal')
        plt.show()
    def showHexBiTmaps(self):
        fig, ax = plt.subplots(constrained_layout=True)
        cc = ax.contourf(self.hexXmap, self.hexYmap, self.TmapG, 100, cmap = 'nipy_spectral', alpha = .5)
        ax.contour(self.hexXmap, self.hexYmap, self.TmapG, 100, cmap = 'nipy_spectral')
        ax.plot(self.path[:,0],self.path[:,1],'r')
#        ax.quiver(self.hexXmap, self.hexYmap,np.cos(self.dirMapG), np.sin(self.dirMapG),scale = 40)
        ax.quiver(self.linkPos[0], self.linkPos[1],np.cos(self.dirLinkG), np.sin(self.dirLinkG),scale = 40)
        cbar = fig.colorbar(cc)
        cbar.set_label('Total Cost To Goal')
        ax.set_xlim([self.xMap[0,0], self.xMap[-1,-1]])
        ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        ax.set_aspect('equal')
        plt.show()  
        fig, ax = plt.subplots(constrained_layout=True)
        cc = ax.contourf(self.hexXmap, self.hexYmap, self.TmapS, 100, cmap = 'nipy_spectral', alpha = .5)
        ax.contour(self.hexXmap, self.hexYmap, self.TmapS, 100, cmap = 'nipy_spectral')
        ax.plot(self.path[:,0],self.path[:,1],'r')
#        ax.quiver(self.hexXmap, self.hexYmap,np.cos(self.dirMapS), np.sin(self.dirMapS),scale = 40)
        ax.quiver(self.linkPos[0], self.linkPos[1],np.cos(self.dirLinkS), np.sin(self.dirLinkS),scale = 40)
        cbar = fig.colorbar(cc)
        cbar.set_label('Total Cost from Start')
        ax.set_xlim([self.xMap[0,0], self.xMap[-1,-1]])
        ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        ax.set_aspect('equal')
        plt.show()   
        fig, ax = plt.subplots(constrained_layout=True)
        cc = ax.contourf(self.hexXmap, self.hexYmap, self.Tmap, 100, cmap = 'nipy_spectral', alpha = .5)
        ax.contour(self.hexXmap, self.hexYmap, self.Tmap, 100, cmap = 'nipy_spectral')
        ax.plot(self.path[:,0],self.path[:,1],'r')
#        ax.quiver(self.hexXmap, self.hexYmap,np.cos(self.dirMapS), np.sin(self.dirMapS),scale = 40)
#        ax.quiver(self.linkPos[0], self.linkPos[1],np.cos(self.dirLinkS), np.sin(self.dirLinkS),scale = 40)
        cbar = fig.colorbar(cc)
        cbar.set_label('Total Cost from Start')
        ax.set_xlim([self.xMap[0,0], self.xMap[-1,-1]])
        ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        ax.set_aspect('equal')
        plt.show()   
    def showResults(self):
        fig, ax = plt.subplots(constrained_layout=True)
        cc = ax.contourf(self.hexXmap, self.hexYmap, self.Tmap, 100, cmap = 'nipy_spectral', alpha = .5)
        ax.contour(self.hexXmap, self.hexYmap, self.Tmap, 100, cmap = 'nipy_spectral')
#        ax.quiver(self.hexXmap, self.hexYmap,np.cos(self.dirMap), np.sin(self.dirMap),scale = 40)
        ax.plot(self.path[:,0],self.path[:,1],'r')
        cbar = fig.colorbar(cc)
        cbar.set_label('Total Cost')
        ax.set_xlim([self.xMap[0,0], self.xMap[-1,-1]])
        ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        ax.set_aspect('equal')
        plt.show()
        
        fig, ax = plt.subplots()
        ax.plot(self.pathTravDist, [x*rad2deg for x in self.pathSlope], linestyle='solid')
        ax.plot(self.pathTravDist, [x*rad2deg for x in self.pathPitch], linestyle='dotted')
        ax.plot(self.pathTravDist, [x*rad2deg for x in self.pathRoll], linestyle='dashed')
        ax.set_aspect('equal')
        ax.legend(('Slope','Pitch','Roll'))
        plt.show()
        
        fig, ax = plt.subplots()
        ax.plot(self.pathTravDist, self.pathComputedTotalCost, linestyle='solid')
        ax.plot(self.pathTravDist, [self.pathEstimatedTotalCost[0]-x for x in self.pathEstimatedTotalCost], linestyle='dotted')
#        ax.set_aspect('equal')
        ax.legend(('Computed', 'Estimated'))
        plt.show()
        
        
    def showPath(self, fig, axes, color, style):
        axes.plot(self.path[:,0],self.path[:,1], color, linestyle=style)
#        axes.scatter(self.linkPos[0], self.linkPos[1], c=color)
        plt.show()
        
    def showPathData(self, opt, fig, axes, color):
        if   opt == 'elevation':
            axes.plot(self.pathTravDist, self.pathElevation, color)
        elif opt == 'slope':
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathSlope], color, linewidth = 2)
        elif opt == 'pitch':
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathPitch], color)
        elif opt == 'roll':
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathRoll], color)
        elif opt == 'cost':
            axes.plot(self.pathTravDist, self.pathCost, color, linewidth = 2)
            axes.set_xlim([self.pathTravDist[0], self.pathTravDist[-1]])
        elif opt == 'total-cost':
            axes.plot(self.pathTravDist, self.pathComputedTotalCost, color, linewidth = 2)
            axes.set_xlim([self.pathTravDist[0], self.pathTravDist[-1]])
        elif opt == 'total-cost-estimated':
            axes.plot(self.pathTravDist, [self.pathEstimatedTotalCost[0]-x for x in self.pathEstimatedTotalCost], color, linewidth = 2, linestyle='dotted')
            axes.set_xlim([self.pathTravDist[0], self.pathTravDist[-1]])
        elif opt == 'total-cost-estimated-filtered':
            axes.plot(self.pathTravDist, [x for x in self.pathEstimatedTotalCost], color, linewidth = 2, linestyle='dotted')
            axes.set_xlim([self.pathTravDist[0], self.pathTravDist[-1]]) 
        elif opt == 'full-orientation':
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathSlope], color, linestyle='solid', linewidth = 2)
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathPitch], color, linestyle='dotted', linewidth = 2)
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathRoll], color, linestyle='dashed', linewidth = 2)
            axes.set_xlim([self.pathTravDist[0], self.pathTravDist[-1]])
        elif opt == 'full-heading':
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathHeading], color, linestyle='solid', linewidth = 2)
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathAspect], color, linestyle='dotted', linewidth = 2)
            axes.plot(self.pathTravDist, [x*rad2deg for x in self.pathBeta], color, linestyle='dashed', linewidth = 2)
            axes.set_xlim([self.pathTravDist[0], self.pathTravDist[-1]])
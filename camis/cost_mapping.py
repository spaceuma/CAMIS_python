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
#============================CAMIS library=====================================
#           Continuous Anisotropic Model for Inclined Surfaces
# -----------------------------------------------------------------------------
#                           cost_mapping.py
#   This file contains a library of python functions dedicated to the 
#   creation of maps with anisotropic cost
#   Notation prefixes:
#     - (nothing)- = referred to input DEM (e.g. elevationMap)
#     - hex- = referred to hexagonal grid (e.g. hexElevationMap) for planning
#     - sq- = referred to square grid (e.g. sqElevationMap) for planning
#==============================================================================

import numpy as np
import math
import sys
try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib import cbook
    from matplotlib import cm
    from matplotlib.colors import LightSource
except:
    raise ImportError('ERROR: matplotlib module could not be imported')

try:
    import scipy.interpolate as interp
    from scipy import signal
    from scipy import ndimage
    import scipy.linalg
except:
    raise ImportError('ERROR: scipy module could not be imported')


import camis.anisotropic_path_planning as ap


try:
    from mayavi import mlab
    isMayavi = True
except:
    isMayavi = False
    print(\
        'WARNING: Mayavi module is not found, some functions cannot be used')
    

from time import time





# Proccessed DEM
class AnisotropicMap:
    def __init__(self, DEM, demRes, planRes, offset, occupancy_radius, 
                 tracking_error):
        # Conversion units
        self.rad2deg = 180/np.pi
        
        # The input Digital Elevation Map in square grid format
        self.elevationMap = DEM
        self.size = DEM.shape
        
        # Field in XYZ of vectors normal to the surface
        self.nx = np.zeros_like(self.elevationMap)
        self.ny = np.zeros_like(self.elevationMap)
        self.nz = np.zeros_like(self.elevationMap)
        
        # The X and Y positions using local frame (local offset is 0,0)
        xMap, yMap = \
          np.meshgrid(np.linspace(0,self.size[1]-1,self.size[1]), \
                      np.linspace(0,self.size[0]-1,self.size[0]))
        self.xMap = xMap*demRes
        self.yMap = yMap*demRes
        self.XYpoints = np.zeros((self.xMap.size,2))
        self.XYpoints[:,0] = self.xMap.flatten()
        self.XYpoints[:,1] = self.yMap.flatten()
        
        # DEM resolution
        self.demRes = demRes
        
        # Resolution of processed map for path planning
        # This is distance between nodes! In hexagonal grid it is not equal to
        # the lateral of the hexagonal cell
        self.planRes = planRes
        
        # If the planning is done in hexagons
        self.hexCellLateral = planRes / math.sqrt(3)
        self.sqCellLateral = math.sqrt(math.sqrt(3)/2*planRes**2)
        
        # This is the global offset from the DEM
        self.offset = offset
        self.xMin = self.xMap[0][0]
        self.yMin = self.yMap[0][0]
        
        # This is the proximity to obstacles, initialized to inf here
        self.proximityMap = np.ones_like(xMap)*np.inf
        
        # Parameters relevant for navigation
        # Radius of the circumference that circumbscripts the robot 2d 
        # projection
        self.occupancy_radius = occupancy_radius
        # Admissible geometric error from the tracking algorithm
        self.tracking_error = tracking_error
        
        init = time()   
        self.computeOccupancyMatrix()
        print('Occupancy Matrix computed in ' + str(time()-init) + ' seconds') 
        
        init = time() 
        self.processElevationMap()
        print('Elevation Map processed in ' + str(time()-init) + ' seconds') 
        
        init = time()
        self.computeObstacleProximityMap()
        print('Obstacle Proximity Map processed in ' + str(time()-init) + 
              ' seconds') 
        
        # Hexagonal navigation maps preparation
        self.computeHexNavMaps()
        
        # Square navigation maps preparation
        self.computeSqNavMaps()


    def computeOccupancyMatrix(self):
        radius = self.occupancy_radius + \
                 self.tracking_error
        print('Occupancy radius is: ', radius, ' meters')
        r = int(radius/self.demRes)
        r = r + 1 - r%2
        y,x = np.ogrid[-r: r+1, -r: r+1]
        print('Occupancy radius in nodes is: ', r, ' meters')
        convMatrix = x**2+y**2 <= r**2
        convMatrix = convMatrix.astype(float)
        self.occupancyMatrix = convMatrix
        self.occupancyMatrixNorm = convMatrix/convMatrix.sum()
        self.radius = radius
    
    def computeObstacleProximityMap(self):
        # We define the forbidden areas
        self.obstacleMap = np.zeros_like(self.elevationMap)
        self.obstacleMap[0,:] = 1
        self.obstacleMap[-1,:] = 1
        self.obstacleMap[:,0] = 1
        self.obstacleMap[:,-1] = 1
        
        self.obstacleMap = ndimage.morphology.binary_dilation(self.obstacleMap,
            structure = self.occupancyMatrix).astype(self.obstacleMap.dtype)
        self.proximityMap = ndimage.morphology.distance_transform_edt(
            1-self.obstacleMap)
        self.proximityMap = self.proximityMap*self.demRes     
        
    def processElevationMap(self):
        processedElevationMap = self.elevationMap
        processedElevationMap = signal.convolve2d(processedElevationMap, \
                                                  self.occupancyMatrixNorm, \
                                                  mode='same', boundary='symm')
        self.processedElevationMap = processedElevationMap
        r = int(self.radius/self.demRes)
        r = np.max((3,r + 1 - r%2))
        r2 = int((r - 1) / 2)
        print('Occupancy radius in nodes is: ', r, ' nodes')
        print('Half occupancy radius in nodes is: ', r2, ' nodes')
        XX,YY = np.meshgrid(np.linspace(-0.5,0.5,r), np.linspace(-0.5,0.5,r))       
        Xarray = XX.flatten()
        Yarray = YY.flatten()
        Zarray = np.zeros_like(Xarray)
        self.aspectX = np.zeros_like(self.elevationMap)
        self.aspectY = np.zeros_like(self.elevationMap)
        self.slope = np.ones_like(self.elevationMap)*np.nan
        
        # Adjusting planes to batches of points
        for j in range(r2,self.processedElevationMap.shape[0] - r2 - 1):
            for i in range(r2,self.processedElevationMap.shape[1] - r2 - 1):
                for l in range(-r2,r2+1):
                    for k in range(-r2,r2+1):
                        Zarray[k+r2 + (l+r2)*r] = \
                            self.processedElevationMap[j+l][i+k] - \
                                self.processedElevationMap[j][i]
                A = np.c_[Xarray, Yarray, np.ones(Xarray.size)]
                C,_,_,_ = scipy.linalg.lstsq(A, Zarray)
                self.nx[j][i] = -C[0] / np.linalg.norm([-C[0], -C[1], 1.0])
                self.ny[j][i] = -C[1] / np.linalg.norm([-C[0], -C[1], 1.0])
                self.nz[j][i] = 1.0 / np.linalg.norm([-C[0], -C[1], 1.0])
                self.aspectX[j][i] = self.nx[j][i] / np.linalg.norm([
                    self.nx[j][i], self.ny[j][i]])
                self.aspectY[j][i] = self.ny[j][i] / np.linalg.norm([
                    self.nx[j][i], self.ny[j][i]])
                self.slope[j][i] = np.abs(np.arccos(self.nz[j][i]))

    def computeHexNavMaps(self):
        init = time()
        self.computeHexGrid()
        print('Hexagonal Grid created in ' + str(time()-init) + ' seconds')
        
        init = time() 
        self.computeHexElevationMap()
        print('Hexagonal Elevation Map created in ' + str(time()-init) 
              + ' seconds') 
        
        init = time() 
        self.computeHexSlopeMaps()
        print('Hexagonal Steepness and Aspect Maps created in ' + 
              str(time()-init) + ' seconds')   
        
        self.hexUsedNodesNumber = np.count_nonzero(~np.isnan(
            self.hexElevationMap))
        self.hexTotalNodesNumber = self.hexElevationMap.size 
        self.hexUsedNodesRatio = self.hexUsedNodesNumber / \
                                 self.hexTotalNodesNumber

    def computeSqNavMaps(self):
        # It would be interesting to detect if demRes = planRes and save up
        # computation
        init = time()
        self.computeSqGrid()
        print('Square Grid created in ' + str(time()-init) + ' seconds')
        
        init = time() 
        self.computeSqElevationMap()
        print('Square Elevation Map created in ' + str(time()-init) 
              + ' seconds') 
        
        init = time() 
        self.computeSqSlopeMaps()
        print('Square Steepness and Aspect Maps created in ' + 
              str(time()-init) + ' seconds')  

        self.sqUsedNodesNumber = np.count_nonzero(~np.isnan(
            self.sqElevationMap))
        self.sqTotalNodesNumber = self.sqElevationMap.size 
        self.sqUsedNodesRatio = self.sqUsedNodesNumber / \
                                self.sqTotalNodesNumber
        
    def computeHexGrid(self):
        DX = self.xMap[-1,-1] - self.xMap[0,0]
        DY = self.yMap[-1,-1] - self.yMap[0,0]
        JMax = math.ceil(2*DY/(math.sqrt(3)*self.planRes))
        IMax = math.ceil((DX+DY/math.sqrt(3))/self.planRes)
        II,JJ = np.meshgrid(np.linspace(0,IMax,IMax+1),
                            np.linspace(0,JMax,JMax+1))
        self.hexXmap = self.xMap[0,0] + self.planRes*(II + .5*JJ) - \
            DY/math.sqrt(3)
        self.hexYmap = self.yMap[0,0] + self.planRes*math.sqrt(3)/2*JJ
        XX,YY = np.meshgrid(np.linspace(0,int(np.ceil(DX)),int(np.ceil(DX))+1),
                            np.linspace(0,int(np.ceil(DY)),int(np.ceil(DY))+1))
        self.xy2J = 2*YY/(np.sqrt(3)*self.planRes)
        self.xy2I = (DY/np.sqrt(3)+ XX)/self.planRes-0.5*self.xy2J  
        
      
    def computeSqGrid(self):
        DX = self.xMap[-1,-1] - self.xMap[0,0]
        DY = self.yMap[-1,-1] - self.yMap[0,0]
        IMax = math.ceil(DX/self.sqCellLateral)
        JMax = math.ceil(DY/self.sqCellLateral)
        II,JJ = np.meshgrid(np.linspace(0,IMax,IMax+1),
                            np.linspace(0,JMax,JMax+1))
        self.sqXmap = self.xMap[0,0] + self.sqCellLateral*II
        self.sqYmap = self.yMap[0,0] + self.sqCellLateral*JJ
        
    def computeHexElevationMap(self):
        self.hexElevationMap = np.ones([self.hexXmap.shape[0],
                                        self.hexXmap.shape[1]])*np.nan        
        self.hexElevationMap[:] = interp.griddata(self.XYpoints,
                                        self.processedElevationMap.flatten(),
                                        (self.hexXmap, self.hexYmap),
                                        method='cubic')
            
    def computeSqElevationMap(self):        
        self.sqElevationMap = interp.griddata(self.XYpoints, \
                                        self.processedElevationMap.flatten(),\
                                        (self.sqXmap, self.sqYmap), \
                                        method='cubic')    

    def computeHexSlopeMaps(self):
        # Hexagonal Steepness Map
        self.hexSlopeMap = np.ones([self.hexElevationMap.shape[0],
                               self.hexElevationMap.shape[1]])*np.nan
        
        self.hexSlopeMap = np.abs(interp.griddata(self.XYpoints,
                                                  self.slope.flatten(),
                                                  (self.hexXmap, self.hexYmap),
                                                  method='linear'))
        
        # Hexagonal Aspect Map
        self.hexAspectMap = np.ones([2,self.hexElevationMap.shape[0],
                                     self.hexElevationMap.shape[1]])*np.nan
        hexAspectMapX = interp.griddata(self.XYpoints, self.aspectX.flatten(),
                                        (self.hexXmap, self.hexYmap),
                                        method='linear')
        hexAspectMapY = interp.griddata(self.XYpoints, self.aspectY.flatten(),
                                        (self.hexXmap, self.hexYmap),
                                        method='linear')
        self.hexAspectMap[0] = hexAspectMapX / np.sqrt(hexAspectMapX**2 + \
                                                       hexAspectMapY**2)
        self.hexAspectMap[1] = hexAspectMapY / np.sqrt(hexAspectMapX**2 + \
                                                       hexAspectMapY**2)
    
    def computeSqSlopeMaps(self):
        # Square Steepness Map
        self.sqSlopeMap = np.ones([self.sqElevationMap.shape[0],
                               self.sqElevationMap.shape[1]])*np.nan
        
        self.sqSlopeMap = np.abs(interp.griddata(self.XYpoints,
                                                  self.slope.flatten(),
                                                  (self.sqXmap, self.sqYmap),
                                                  method='linear'))
        
        # Square Aspect Map
        self.sqAspectMap = np.ones([2,self.sqElevationMap.shape[0],
                                     self.sqElevationMap.shape[1]])*np.nan
        sqAspectMapX = interp.griddata(self.XYpoints, self.aspectX.flatten(),
                                        (self.sqXmap, self.sqYmap),
                                        method='linear')
        sqAspectMapY = interp.griddata(self.XYpoints, self.aspectY.flatten(),
                                        (self.sqXmap, self.sqYmap),
                                        method='linear')
        self.sqAspectMap[0] = sqAspectMapX / np.sqrt(sqAspectMapX**2 + \
                                                       sqAspectMapY**2)
        self.sqAspectMap[1] = sqAspectMapY / np.sqrt(sqAspectMapX**2 + \
                                                       sqAspectMapY**2)
    
    def computeVecCostMap(self, costModel):
        
        # Imported the CAMIS Cost Model
        self.costModel = costModel
        Cmax = self.costModel.limit_cost
        
        init = time()        
        hexVectorialData = costModel.getVectorialCostMap(
            self.rad2deg*self.hexSlopeMap)
        print('Elapsed time to compute the Vectorial Data: '+str(time()-init))
        
        init = time()        
        sqVectorialData = costModel.getVectorialCostMap(
            self.rad2deg*self.sqSlopeMap)
        print('Elapsed time to compute the Vectorial Data: '+str(time()-init)) 
        
        
        # Creating the obstacle proximity maps
        init = time() 
        self.computeHexObstacleProximityMap()
        
        hexObstacleMask = self.hexProximityMap <= \
            self.costModel.occupancy_radius + \
            self.costModel.tracking_error + sys.float_info.epsilon
            
        self.hexAnisotropyMap = hexVectorialData[0][:][:]
        self.hexAnisotropyMap[hexObstacleMask] = np.inf
        
        self.hexVCMap = np.zeros([4,self.hexAnisotropyMap.shape[0],
                                  self.hexAnisotropyMap.shape[1]])
        self.hexVCMap[0] = hexVectorialData[1][:][:]
        self.hexVCMap[0][np.where(self.hexProximityMap[:]<self.radius)] = \
                                                                        Cmax**3
        self.hexVCMap[0][hexObstacleMask] = np.inf
        self.hexVCMap[1] = hexVectorialData[2][:][:]
        self.hexVCMap[1][np.where(self.hexProximityMap[:]<self.radius)] = \
                                                                        Cmax**3
        self.hexVCMap[1][hexObstacleMask] = np.inf
        self.hexVCMap[2] = hexVectorialData[3][:][:]
        self.hexVCMap[2][np.where(self.hexProximityMap[:]<self.radius)] = 0
        self.hexVCMap[2][hexObstacleMask] = 0
        self.hexVCMap[3] = hexVectorialData[4][:][:]
        self.hexVCMap[3][np.where(self.hexProximityMap[:]<self.radius)] = 0
        self.hexVCMap[3][hexObstacleMask] = 0
        print('The Hexagonal Cost Map is created in ' + str(time()-init)  + 
              ' seconds')
        
        
        init = time() 
        self.computeSqObstacleProximityMap()

        sqObstacleMask = self.sqProximityMap <= \
            self.costModel.occupancy_radius + \
            self.costModel.tracking_error + sys.float_info.epsilon
            
        self.sqAnisotropyMap = sqVectorialData[0][:][:]
        self.sqAnisotropyMap[sqObstacleMask] = np.inf
        
        self.sqVCMap = np.zeros([4,self.sqAnisotropyMap.shape[0],
                                  self.sqAnisotropyMap.shape[1]])
        self.sqVCMap[0] = sqVectorialData[1][:][:]
        self.sqVCMap[0][np.where(self.sqProximityMap[:]<self.radius)] = \
                                                                        Cmax**3
        self.sqVCMap[0][sqObstacleMask] = np.inf
        self.sqVCMap[1] = sqVectorialData[2][:][:]
        self.sqVCMap[1][np.where(self.sqProximityMap[:]<self.radius)] = \
                                                                        Cmax**3
        self.sqVCMap[1][sqObstacleMask] = np.inf
        self.sqVCMap[2] = sqVectorialData[3][:][:]
        self.sqVCMap[2][np.where(self.sqProximityMap[:]<self.radius)] = 0
        self.sqVCMap[2][sqObstacleMask] = 0
        self.sqVCMap[3] = sqVectorialData[4][:][:]
        self.sqVCMap[3][np.where(self.sqProximityMap[:]<self.radius)] = 0
        self.sqVCMap[3][sqObstacleMask] = 0
        print('The Square Cost Map is created in ' + str(time()-init)  + 
              ' seconds')
        
        
    def computeHexObstacleProximityMap(self):
        self.hexProximityMap = interp.griddata(self.XYpoints, 
                                               self.proximityMap.flatten(),\
                                               (self.hexXmap, self.hexYmap), 
                                               method='nearest')
        self.hexProximityMap[np.where(np.isnan(self.hexProximityMap))] = 0.0
        
        
    def computeSqObstacleProximityMap(self):
        self.sqProximityMap = interp.griddata(self.XYpoints, 
                                               self.proximityMap.flatten(),\
                                               (self.sqXmap, self.sqYmap), 
                                               method='nearest')
        self.sqProximityMap[np.where(np.isnan(self.sqProximityMap))] = 0.0   
        
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
        ap.computeTmap(self.hexVCMap, self.hexAspectMap, self.hexAnisotropyMap,
                       ijGoal, ijStart, self.hexXmap, self.hexYmap, 
                       self.planRes)
        elapsedTime = time()-init
        print('Elapsed time to compute the Total Cost Map: '+str(elapsedTime))
        
        
#        startWaypoint = IJ2XY[:,start[1],start[0]]
#        goalWaypoint = IJ2XY[:,goal[1],goal[0]]
        self.IJ2XY = IJ2XY
        self.XY2IJ = XY2IJ
        self.dirMap = dirMap
        path,uu = ap.getPath(dirMap, IJ2XY, XY2IJ, start, goal, self.xMin, self.yMin, self.planRes)
        self.path = np.asarray(path)
        self.getPathData()
        return elapsedTime
    
    # Executing BiOUM to compute the path faster
    def executeSqBiPlanning(self, goal, start, nbUpdate, anisoSearch, 
                            dirPolicy):
        init = time()
        ijStart = np.round([(start[0] - 
                            self.sqXmap[0][0])/self.sqCellLateral,(start[1] - 
                            self.sqYmap[0][0])/self.sqCellLateral]).astype(int) 
        ijGoal = np.round([(goal[0] - 
                            self.sqXmap[0][0])/self.sqCellLateral,(goal[1] - 
                            self.sqYmap[0][0])/self.sqCellLateral]).astype(int) 
        
        print('The start indexes are ' + str(ijStart) + 
              ' and the position is ' + str(self.sqXmap[ijStart[1]][ijStart[0]]) + 
              ','  + str(self.sqYmap[ijStart[1]][ijStart[0]]))
        
        self.sqTmapG, self.sqTmapS, self.sqDirMapG, self.sqDirMapS, \
            self.sqLinkNode, stateMapG, stateMapS, self.dirLinkG, \
                self.dqDirLinkS = \
        ap.computeBiTmap(self.sqVCMap, self.sqAspectMap, 
                         self.sqAnisotropyMap, ijGoal, ijStart, self.sqXmap, 
                         self.sqYmap, self.sqCellLateral, gridtype='sq',
                         nbUpdate = nbUpdate, anisoSearch = anisoSearch,
                         dirPolicy = dirPolicy)
        
        print('Elapsed time to compute the Total Cost Map: ' + 
              str(time()-init)  +  ' seconds')
        if dirPolicy == 'bidir':
            TlinkG = self.sqTmapG[self.sqLinkNode[1],self.sqLinkNode[0]]
            TlinkS = self.sqTmapS[self.sqLinkNode[1],self.sqLinkNode[0]]
            self.sqTmap = np.zeros_like(self.sqTmapG)
            self.sqTmap[:] = self.sqTmapG
            self.sqTmap[np.where(self.sqTmapS != np.inf)] = TlinkS + TlinkG - self.sqTmapS[np.where(self.sqTmapS != np.inf)]
            self.sqTmap[np.where(self.sqTmap < 0.0)] = np.inf
            
            self.sqLinkPos = np.empty([2],dtype=float)
            self.sqLinkPos[:] = [self.sqXmap[self.sqLinkNode[1],
                                             self.sqLinkNode[0]],
                                 self.sqYmap[self.sqLinkNode[1],
                                             self.sqLinkNode[0]]]
            
            self.sqPathG,uu = ap.getSqPath(self.sqDirMapG, self.sqLinkPos, goal, 
                                         self.xMin, self.yMin, self.sqCellLateral)
            self.sqPathS,uu = ap.getSqPath(self.sqDirMapS, self.sqLinkPos, start, 
                                         self.xMin, self.yMin, self.sqCellLateral)
            self.sqPath = np.concatenate((np.flipud(np.asarray(self.sqPathS)),
                                        np.asarray(self.sqPathG)))
        if dirPolicy == 'goal':
            self.sqTmap = self.sqTmapG
            path,uu = ap.getSqPath(self.sqDirMapG, start, goal, 
                                         self.xMin, self.yMin, 
                                         self.sqCellLateral)
            self.sqPath = np.asarray(path)
        # self.getSqPathData()
    # nbUpdate = True, False
    # anisoSearch = 'single', 'double'
    # dirPolicy = 'start','goal','bidir'
    def executeHexBiPlanning(self, goal, start, nbUpdate, anisoSearch, 
                             dirPolicy):
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
        self.TmapG, self.TmapS, self.dirMapG, self.dirMapS, nodeLink,\
        stateMapG, stateMapS, self.dirLinkG, self.dirLinkS, self.numUpdates = \
        ap.computeBiTmap(self.hexVCMap, self.hexAspectMap, 
                         self.hexAnisotropyMap, ijGoal, ijStart, self.hexXmap, 
                         self.hexYmap, self.planRes, gridtype='hex',
                         nbUpdate = nbUpdate, anisoSearch = anisoSearch,
                         dirPolicy = dirPolicy)
        elapsedTime = time()-init
        print('Elapsed time to compute the Total Cost Map: '+str(elapsedTime))
        self.linkNode = nodeLink
        self.linkPos = IJ2XY[:,nodeLink[1],nodeLink[0]]
        self.IJ2XY = IJ2XY
        self.XY2IJ = XY2IJ
        if dirPolicy == 'bidir':
            self.pathG,uu = ap.getPath(self.dirMapG, IJ2XY, XY2IJ, self.linkPos, 
                                       goal, self.xMin, self.yMin, self.planRes)
            self.pathS,uu = ap.getPath(self.dirMapS, IJ2XY, XY2IJ, self.linkPos, 
                                       start, self.xMin, self.yMin, self.planRes)
        
            TlinkG = self.TmapG[nodeLink[1],nodeLink[0]]
            TlinkS = self.TmapS[nodeLink[1],nodeLink[0]]
        
            self.pathS = self.pathS[1:-1]
            
            self.path = np.concatenate((np.flipud(np.asarray(self.pathS)),
                                        np.asarray(self.pathG)))
            self.Tmap = np.zeros_like(self.TmapG)
            self.Tmap[:] = self.TmapG
            self.Tmap[np.where(self.TmapS != np.inf)] = TlinkS + TlinkG - self.TmapS[np.where(self.TmapS != np.inf)]
            self.Tmap[np.where(self.Tmap < 0.0)] = np.inf
        self.getPathData()
        self.elapsedTime = elapsedTime
        return elapsedTime
        
    def getBeta(self, xPos, yPos, heading):
        beta = []
        for index, waypoint in enumerate(xPos):
            try:
                aX = ap.interpolatePoint([xPos[index]/self.demRes,yPos[index]],self.aspectX)
                aY = ap.interpolatePoint([xPos[index]/self.demRes,yPos[index]],self.aspectY)
            except:
                print('ERROR at index ' + str(index))
                print('X = ' + str(xPos[index]))
                print('Y = ' + str(yPos[index]))
                aX = np.nan
                aY = np.nan
            aspect = np.arctan2(aY,aX)
            b = np.arccos(np.cos(heading[index])*np.cos(aspect)+np.sin(heading[index])*np.sin(aspect))
            crossDirection = np.sin(heading[index])*np.cos(aspect)-np.cos(heading[index])*np.sin(aspect)
            if crossDirection >= 0.:
                beta.append(b)
            else:
                beta.append(-b)
        return beta
    
    def getSlope(self, xPos, yPos, heading):
        slope = []
        for index, waypoint in enumerate(xPos):
            try:
                s = ap.interpolatePoint([xPos[index]/self.demRes,yPos[index]/self.demRes],self.slope)
            except:
                print('ERROR at index ' + str(index))
                print('X = ' + str(xPos[index]))
                print('Y = ' + str(yPos[index]))
                s = np.nan
            slope.append(s)
        return slope
       
    def getSqPathData(self):
        pathElevation = []
        for index, waypoint in enumerate(self.sqPath):
            pathElevation.append(ap.interpolatePoint(
                waypoint/self.sqCellLateral,self.elevationMap))
        self.pathElevation = pathElevation
        
    def getPathData(self):
        pathElevation = []
        pathSlope = []
        pathAspectX = []
        pathAspectY = []
        pathTravDist = []
        pathPitch = []
        pathHeading = []
        pathAspect = []
        pathCost = []
        pathCostwithRisk = []
        pathSegment = []
        pathEstimatedTotalCost = []
        pathComputedTotalCost = []
        pathComputedTotalCostwithRisk = []
        for index, waypoint in enumerate(self.path):
            pathElevation.append(ap.getTriInterpolation(waypoint,self.hexElevationMap,self.XY2IJ,self.IJ2XY,self.planRes, self.xMin, self.yMin))
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
                A = self.path[index] - self.path[index-1]
                pathHeading.append(np.arctan2(A[1],A[0]))
                pathSegment.append(np.linalg.norm(A)/2)
            elif index == 0:
                A = self.path[index+1] - self.path[index]
                pathSegment.append(np.linalg.norm(A)/2)
                pathHeading.append(np.arctan2(A[1],A[0]))
            else:
                A1 = self.path[index+1] - self.path[index]
                A2 = self.path[index] - self.path[index-1]
                pathSegment.append((np.linalg.norm(A1) + \
                                    np.linalg.norm(A2))/2)
                pathHeading.append(np.arctan2(A1[1]+A2[1],A1[0]+A2[0]))
            pathCost.append(self.costModel.getRawCost(self.rad2deg*pathSlope[index],pathAspect[index],pathHeading[index]))
            pathCostwithRisk.append(self.costModel.getCost(self.rad2deg*pathSlope[index],pathAspect[index],pathHeading[index]))
            pathEstimatedTotalCost.append(ap.getTriLowest(waypoint,self.Tmap,self.XY2IJ,self.IJ2XY,self.planRes, self.xMin, self.yMin))
            if index == 0:
                pathComputedTotalCost.append(0)
                pathComputedTotalCostwithRisk.append(0)
            else:
                pathComputedTotalCost.append(pathComputedTotalCost[index-1] + pathCost[index]*pathSegment[index])
                pathComputedTotalCostwithRisk.append(pathComputedTotalCostwithRisk[index-1] + pathCostwithRisk[index]*pathSegment[index])
#                pathComputedTotalCost.append(pathComputedTotalCost[index-1]+0.5*(pathCost[index-1]+pathCost[index])*pathSegment[index-1])
            
        b = np.arccos(np.cos(pathHeading)*np.cos(pathAspect)+np.sin(pathHeading)*np.sin(pathAspect))
        crossDirection = np.sin(pathHeading)*np.cos(pathAspect)-np.cos(pathHeading)*np.sin(pathAspect)
        for i,x in enumerate(b):
            if crossDirection[i] < 0.0:
                b[i] = - b[i]
        
        pathPitch = np.arccos(np.cos(pathSlope)/np.sqrt(np.cos(b)**2+np.cos(pathSlope)**2*np.sin(b)**2))
        for i,x in enumerate(pathPitch):
            if b[i] > np.pi/2 or b[i] < -np.pi/2:
                pathPitch[i] = - pathPitch[i]
        self.pathPitch = pathPitch
        self.pathRoll = np.arccos(np.cos(pathSlope)/np.cos(pathPitch))
        self.pathElevation = pathElevation
        self.pathSlope = pathSlope
        self.pathTravDist = pathTravDist
        self.pathAspect = pathAspect
        self.pathHeading = pathHeading
        self.pathBeta = b
        self.pathCost = pathCost
        self.pathCostwithRisk = pathCostwithRisk
        self.pathSegment = pathSegment
        self.pathComputedTotalCost = pathComputedTotalCost
        self.pathComputedTotalCostwithRisk = pathComputedTotalCostwithRisk
        self.pathEstimatedTotalCost = pathEstimatedTotalCost
        
    def show3dDEM(self):
        if isMayavi:
            mlab.figure(size=(800, 400),bgcolor=(1,1,1), fgcolor=(0.,0.,0.))
            mlab.mesh(np.rot90(self.xMap), np.rot90(self.yMap), np.rot90(self.processedElevationMap), colormap='gist_earth')
            axes = mlab.axes(x_axis_visibility = True, y_axis_visibility = True, color = (0,0,0), 
                      xlabel = 'X-axis [m]', ylabel = 'Y-axis [m]', zlabel = 'Z-axis [m]', extent = (0,25,0,80,53,57))
            axes.label_text_property.font_family = 'times'
            axes.label_text_property.font_size = 10
            axes.axes.label_format='%.1f'
            axes.axes.axis_label_text_property.bold = 0
#            mlab.view(-59, 58, 1773, [-.5, -.5, 512])
        else:
            raise ImportError('show3dDEM: Mayavi is not available')
            
    def showMap(self, opt, fig, axes):
        if   opt == 'elevation':
            cc = axes.contourf(self.xMap, self.yMap, self.elevationMap, 50, cmap = cm.gist_earth, extend='both')
            axes.contour(self.xMap, self.yMap, self.elevationMap, 10, colors = 'k', alpha=.3)
#            cbar = fig.colorbar(cc, orientation="horizontal",fraction=0.046, pad=0.04)
            cbar = fig.colorbar(cc)
            cbar.set_label('Elevation (m)')
        elif opt == 'conv-slope':
#            cc = axes.contourf(self.xMap, self.yMap, 180/3.1416*np.arccos(self.nz), 50, cmap="nipy_spectral")
            cc = axes.contourf(self.xMap, self.yMap, self.aspectY, 50, cmap = cm.gist_earth, extend='both')
#            axes.contour(self.xMap, self.yMap, 180/3.1416*np.arccos(self.nz), 10, colors = 'k', alpha=.3)
#            cc.set_clim(0,30.0)
#            cbar = fig.colorbar(cc, orientation="horizontal",fraction=0.046, pad=0.04)
            cbar = fig.colorbar(cc)
            cbar.set_label('Slope')
        elif opt == 'hex-elevation':
            cc = axes.contourf(self.hexXmap, self.hexYmap, self.hexElevationMap, 50, cmap = cm.gist_earth, extend='both')
            axes.contour(self.hexXmap, self.hexYmap, self.rad2deg*self.hexElevationMap, 20, colors = 'k', alpha=.3)
            cbar = fig.colorbar(cc)
            cbar.set_label('Elevation (m)')
            axes.set_xlim([self.xMap[0,0], self.xMap[-1,-1]])
            axes.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
            axes.set_aspect('equal')
            axes.set_xlabel('X-axis [m]')
            axes.set_ylabel('Y-axis [m]')
        elif opt == 'hex-slope':
            cc = axes.scatter(self.hexXmap, self.hexYmap, c = self.rad2deg*self.hexSlopeMap, cmap="nipy_spectral",s=20)
            axes.set_aspect('equal')
            cbar = fig.colorbar(cc)
            cbar.set_label('Slope (deg)')
            axes.set_aspect('equal')
            axes.set_xlim([self.xMap[0,2], self.xMap[-1,-4]])
            axes.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
            axes.set_xlabel('X-axis [m]')
            axes.set_ylabel('Y-axis [m]')
            plt.show()
        elif opt == 'proximity':
            cc = axes.contourf(self.xMap, self.yMap, self.proximityMap, 100, cmap = 'nipy_spectral')
            cc.set_clim(0,5.0)
            cbar = fig.colorbar(cc)
            cbar.set_label('Proximity to obstacles (m)')
        else:
            raise ValueError(opt + ' is not a valid option')
        axes.set_aspect('equal')
        plt.show()
#    @classmethod
#    def from_pdem(cls, parentPdem, costModel):
#        return cls(parentPdem.elevationMap, parentPdem.demRes,\
#                   parentPdem.planRes, parentPdem.offset, costModel)

    def showHexVecCostMap(self, index):
        fig, ax = plt.subplots()
        chosenVCmap = self.hexVCMap[index]
        chosenVCmap[np.where(chosenVCmap == np.inf)] = np.nan 
        cc = ax.contourf(self.hexXmap, self.hexYmap, chosenVCmap, 100, cmap = 'nipy_spectral')
#        cc.set_clim(0,200.0)
        cbar = fig.colorbar(cc)
        ax.set_aspect('equal')
        plt.show()
        
    def showHexAnisotropyMap(self):
        fig, ax = plt.subplots()
        cc = ax.contourf(self.hexXmap, self.hexYmap, self.hexAnisotropyMap, 100, cmap = 'nipy_spectral')
        cbar = fig.colorbar(cc)
        cbar.set_label('Anisotropy')
        ax.set_aspect('equal')
        plt.show()
       
    def showElevationMap(self,fontsize):
        fig, ax = plt.subplots(constrained_layout=True)
        cc = ax.contourf(self.xMap, self.yMap, self.elevationMap, 
                         50, cmap = cm.gist_earth, extend='both')
        ax.contour(self.xMap, self.yMap, self.hexElevationMap, 
                   20, colors = 'k', alpha=.3)
        # cc = ax.scatter(self.hexXmap, self.hexYmap, 
        #                 c = self.hexElevationMap, cmap=cm.gist_earth,s=20)
        ax.set_xlim([self.xMap[0,2], self.xMap[-1,-4]])
        ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        cbar = fig.colorbar(cc)
        cbar.set_label('Elevation [m]', fontsize = fontsize)
        ax.set_aspect('equal')
        
        ax.set_xlabel('X-axis [m]', fontsize = fontsize)
        ax.set_ylabel('Y-axis [m]', fontsize = fontsize)
        
    def showHexElevationMap(self, fontsize):
        self.showNavElevationMap('hex', fontsize)
        
    def showSqElevationMap(self, fontsize):
        self.showNavElevationMap('sq', fontsize)
    
    def showNavElevationMap(self, option, fontsize):
        fig, ax = plt.subplots(constrained_layout=True)
        if option == 'sq':
            cc = ax.scatter(self.sqXmap, self.sqYmap, c = self.sqElevationMap, 
                            cmap=cm.gist_earth,s=20)
            ax.set_xlim([self.xMap[0,0], self.xMap[-1,-1]])
            ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        elif option == 'hex':
            cc = ax.contourf(self.hexXmap, self.hexYmap, self.hexElevationMap, 
                             50, cmap = cm.gist_earth, extend='both')
            ax.contour(self.hexXmap, self.hexYmap, self.hexElevationMap, 
                       10, colors = 'k', alpha=.3)
            # cc = ax.scatter(self.hexXmap, self.hexYmap, 
            #                 c = self.hexElevationMap, cmap=cm.gist_earth,s=20)
            ax.set_xlim([self.xMap[0,2], self.xMap[-1,-4]])
            ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        cbar = fig.colorbar(cc)
        cbar.set_label('Elevation [m]', fontsize = fontsize)
        ax.set_aspect('equal')
        
        ax.set_xlabel('X-axis [m]', fontsize = fontsize)
        ax.set_ylabel('Y-axis [m]', fontsize = fontsize)

    def showSqSlopeGradientMap(self, fontsize, maxGradient = 0.0):
        fig, ax = plt.subplots(constrained_layout=True)
        cc = ax.scatter(self.sqXmap, self.sqYmap, c = 
                        self.rad2deg*self.sqSlopeMap, 
                        cmap="nipy_spectral",s=20)
        if maxGradient > 0.0:
            cc.set_clim(0,maxGradient)
        ax.set_aspect('equal')
        cbar = fig.colorbar(cc)
        cbar.set_label('Slope Gradient [deg]', fontsize = fontsize)
        ax.set_aspect('equal')
        ax.set_xlim([self.xMap[0,0], self.xMap[-1,-1]])
        ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        ax.set_xlabel('X-axis [m]', fontsize = fontsize)
        ax.set_ylabel('Y-axis [m]', fontsize = fontsize)
        plt.show()    
    
    def showHexSteepnessMap(self, fontsize):
        fig, ax = plt.subplots(constrained_layout=True)
        cc = ax.scatter(self.hexXmap, self.hexYmap, c = self.rad2deg*self.hexSlopeMap, cmap="nipy_spectral",s=20)
        ax.set_aspect('equal')
        cbar = fig.colorbar(cc)
        cbar.set_label('Slope (deg)', fontsize = fontsize)
        ax.set_aspect('equal')
        ax.set_xlim([self.xMap[0,2], self.xMap[-1,-4]])
        ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        ax.set_xlabel('X-axis [m]', fontsize = fontsize)
        ax.set_ylabel('Y-axis [m]', fontsize = fontsize)
        plt.show()
        
        
    def showHexAspectMap(self):
        fig, ax = plt.subplots(constrained_layout=True)
        cc = ax.scatter(self.hexXmap, self.hexYmap, c = self.rad2deg*np.arctan2(self.hexAspectMap[1],self.hexAspectMap[0]), cmap="hsv",s=20)
        ax.set_aspect('equal')
        
#        fig, ax = plt.subplots(constrained_layout=True)
#        cc = ax.contourf(self.hexXmap, self.hexYmap, self.rad2deg*self.hexSlopeMap, levels = levels, cmap = 'nipy_spectral', extend = 'max')
#        ax.quiver(self.hexXmap, self.hexYmap,self.hexAspectMap[0], self.hexAspectMap[1],scale = 100)
        cbar = fig.colorbar(cc)
        cbar.set_label('Aspect (deg)')
        ax.set_aspect('equal')
        ax.set_xlim([self.xMap[0,2], self.xMap[-1,-4]])
        ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        ax.set_xlabel('X-axis [m]')
        ax.set_ylabel('Y-axis [m]')
        plt.show()
        
    def showSqAspectMap(self):
        fig, ax = plt.subplots(constrained_layout=True)
        cc = ax.scatter(self.sqXmap, self.sqYmap, 
                        c = self.rad2deg*np.arctan2(self.sqAspectMap[1],
                                                    self.sqAspectMap[0]), 
                        cmap="hsv",s=20)
        ax.set_aspect('equal')
        
        cbar = fig.colorbar(cc)
        cbar.set_label('Aspect (deg)')
        ax.set_aspect('equal')
        ax.set_xlim([self.xMap[0,0], self.xMap[-1,-1]])
        ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        ax.set_xlabel('X-axis [m]')
        ax.set_ylabel('Y-axis [m]')
        plt.show()
        
    def showHexBiTmaps(self):
#         fig, ax = plt.subplots(constrained_layout=True)
#         cc = ax.contourf(self.hexXmap, self.hexYmap, self.TmapG, 100, cmap = 'nipy_spectral', alpha = .5)
#         ax.contour(self.hexXmap, self.hexYmap, self.TmapG, 100, cmap = 'nipy_spectral')
#         ax.quiver(self.hexXmap, self.hexYmap,np.cos(self.dirMapG), np.sin(self.dirMapG),scale = 40)
#         ax.plot(self.path[:,0],self.path[:,1],'r')
# #        ax.quiver(self.hexXmap, self.hexYmap,np.cos(self.dirMapG), np.sin(self.dirMapG),scale = 40)
# #        ax.quiver(self.linkPos[0], self.linkPos[1],np.cos(self.dirLinkG), np.sin(self.dirLinkG),scale = 40)
#         cbar = fig.colorbar(cc)
#         cbar.set_label('Total Cost To Goal')
#         ax.set_xlim([self.xMap[0,0], self.xMap[-1,-1]])
#         ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
#         ax.set_aspect('equal')
#         plt.show()  
#         fig, ax = plt.subplots(constrained_layout=True)
#         cc = ax.contourf(self.hexXmap, self.hexYmap, self.TmapS, 100, cmap = 'nipy_spectral', alpha = .5)
#         ax.contour(self.hexXmap, self.hexYmap, self.TmapS, 100, cmap = 'nipy_spectral')
#         ax.plot(self.path[:,0],self.path[:,1],'r')
# #        ax.quiver(self.hexXmap, self.hexYmap,np.cos(self.dirMapS), np.sin(self.dirMapS),scale = 40)
#         ax.quiver(self.linkPos[0], self.linkPos[1],np.cos(self.dirLinkS), np.sin(self.dirLinkS),scale = 40)
#         cbar = fig.colorbar(cc)
#         cbar.set_label('Total Cost from Start')
#         ax.set_xlim([self.xMap[0,0], self.xMap[-1,-1]])
#         ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
#         ax.set_aspect('equal')
#         plt.show()   
        fig, ax = plt.subplots(figsize=(5.5,4.5),constrained_layout=True)
        cc = ax.contourf(self.hexXmap, self.hexYmap, self.Tmap, 100, cmap = 'nipy_spectral', alpha = .5)
        ax.contour(self.hexXmap, self.hexYmap, self.Tmap, 100, cmap = 'nipy_spectral')
        ax.plot(self.path[:,0],self.path[:,1],'r')
        # ax.quiver(self.hexXmap, self.hexYmap,np.cos(self.dirMapS), np.sin(self.dirMapS),scale = 40)
        # ax.quiver(self.hexXmap, self.hexYmap,np.cos(self.dirMapG), np.sin(self.dirMapG),scale = 40)
#        ax.quiver(self.linkPos[0], self.linkPos[1],np.cos(self.dirLinkS), np.sin(self.dirLinkS),scale = 40)
        cbar = fig.colorbar(cc)
        cbar.set_label('T($x_{ij},x_o$) [As/m]', fontsize = 14)
        ax.set_xlim([self.xMap[0,0], self.xMap[-1,-1]])
        ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        ax.set_ylabel('Y-axis [m]', fontsize = 14)
        ax.set_xlabel('X-axis [m]', fontsize = 14)
        ax.set_aspect('equal')
        plt.show()   
    
    
    def showSqBiTmaps(self):
        fig, ax = plt.subplots(constrained_layout=True)
        cc = ax.contourf(self.sqXmap, self.sqYmap, self.sqTmap, 100, 
                         cmap = 'nipy_spectral', alpha = .5)
        ax.contour(self.sqXmap, self.sqYmap, self.sqTmap, 100, 
                   cmap = 'nipy_spectral')
        ax.plot(self.sqPath[:,0],self.sqPath[:,1],'r')
        cbar = fig.colorbar(cc)
        cbar.set_label('Total Cost')
        ax.set_xlim([self.xMap[0,0], self.xMap[-1,-1]])
        ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        ax.set_aspect('equal')
        plt.show()
    
    def showSqBiHeading(self):
        fig, ax = plt.subplots(constrained_layout=True)
        cc = ax.scatter(self.sqXmap, self.sqYmap, 
                        c = self.rad2deg*self.sqDirMapS, 
                        cmap="hsv",s=20)
        ax.set_aspect('equal')
        
        cbar = fig.colorbar(cc)
        cbar.set_label('Heading (deg)')
        ax.set_aspect('equal')
        ax.set_xlim([self.xMap[0,0], self.xMap[-1,-1]])
        ax.set_ylim([self.yMap[0,0], self.yMap[-1,-1]])
        ax.set_xlabel('X-axis [m]')
        ax.set_ylabel('Y-axis [m]')
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
        ax.plot(self.pathTravDist, np.asarray(self.pathSlope)*self.rad2deg, linestyle='solid')
        ax.plot(self.pathTravDist, np.asarray(self.pathPitch)*self.rad2deg, linestyle='dotted')
        ax.plot(self.pathTravDist, np.asarray(self.pathRoll)*self.rad2deg, linestyle='dashed')
        ax.set_aspect('equal')
        ax.legend(('Slope','Pitch','Roll'))
        plt.show()
        
        fig, ax = plt.subplots()
        ax.plot(self.pathTravDist, self.pathComputedTotalCost, linestyle='solid')
        ax.plot(self.pathTravDist, [self.pathEstimatedTotalCost[0]-x for x in self.pathEstimatedTotalCost], linestyle='dotted')
#        ax.set_aspect('equal')
        ax.legend(('Computed', 'Estimated'))
        plt.show()
        
#        fig, ax = plt.subplots()
#        xx = self.hexXmap.flatten()
#        yy = self.hexYmap.flatten()
#        tt = self.TmapG.flatten()
#        cc = np.cos(self.dirMapG).flatten()
#        ss = np.sin(self.dirMapG).flatten()
#        strm = ax.streamplot(xx, yy, cc, ss, linewidth=2, cmap='autumn')
#        fig.colorbar(strm.lines)
#        ax.set_title('T Map G')
#        
        
    def showPath(self, fig, axes, color, style, alpha=1.0):
        axes.plot(self.path[:,0],self.path[:,1], color, linestyle=style, alpha = alpha)
#        axes.scatter(self.linkPos[0], self.linkPos[1], c=color)
        plt.show()
        
    def showPathData(self, opt, fig, axes, color, style):
        if   opt == 'elevation':
            axes.plot(self.pathTravDist, self.pathElevation, color)
        elif opt == 'segment':
            axes.plot(self.pathTravDist, self.pathSegment, color, linestyle=style)
        elif opt == 'slope':
            axes.plot(self.pathTravDist, [x*self.rad2deg for x in self.pathSlope], color, linewidth = 2)
        elif opt == 'heading':
            axes.plot(self.pathTravDist, [x*self.rad2deg for x in self.pathHeading], color, linestyle=style)
        elif opt == 'beta':
            axes.plot(self.pathTravDist, [x*self.rad2deg for x in self.pathBeta], color, linestyle=style)
        elif opt == 'pitch':
            axes.plot(self.pathTravDist, [x*self.rad2deg for x in self.pathPitch], color)
        elif opt == 'roll':
            axes.plot(self.pathTravDist, [x*self.rad2deg for x in self.pathRoll], color)
        elif opt == 'cost':
            axes.plot(self.pathTravDist, self.pathCost, color, linewidth = 2)
            axes.set_xlim([self.pathTravDist[0], self.pathTravDist[-1]])
        elif opt == 'total-cost':
            axes.plot(self.pathTravDist, self.pathComputedTotalCost, color, linestyle=style, linewidth = 2)
            axes.set_xlim([self.pathTravDist[0], self.pathTravDist[-1]])
        elif opt == 'total-cost-estimated':
            axes.plot(self.pathTravDist, [self.pathEstimatedTotalCost[0]-x for x in self.pathEstimatedTotalCost], color, linewidth = 2, linestyle='dotted')
            axes.set_xlim([self.pathTravDist[0], self.pathTravDist[-1]])
        elif opt == 'total-cost-estimated-filtered':
            axes.plot(self.pathTravDist, [x for x in self.pathEstimatedTotalCost], color, linewidth = 2, linestyle='dotted')
            axes.set_xlim([self.pathTravDist[0], self.pathTravDist[-1]]) 
        elif opt == 'full-orientation':
            axes.plot(self.pathTravDist, [x*self.rad2deg for x in self.pathSlope], color, linestyle='solid', linewidth = 2)
            axes.plot(self.pathTravDist, [x*self.rad2deg for x in self.pathPitch], color, linestyle='dotted', linewidth = 2)
            axes.plot(self.pathTravDist, [x*self.rad2deg for x in self.pathRoll], color, linestyle='dashed', linewidth = 2)
            axes.set_xlim([self.pathTravDist[0], self.pathTravDist[-1]])
        elif opt == 'full-heading':
            axes.plot(self.pathTravDist, [x*self.rad2deg for x in self.pathHeading], color, linestyle='solid', linewidth = 2)
            axes.plot(self.pathTravDist, [x*self.rad2deg for x in self.pathAspect], color, linestyle='dotted', linewidth = 2)
            axes.plot(self.pathTravDist, [x*self.rad2deg for x in self.pathBeta], color, linestyle='dashed', linewidth = 2)
            axes.set_xlim([self.pathTravDist[0], self.pathTravDist[-1]])
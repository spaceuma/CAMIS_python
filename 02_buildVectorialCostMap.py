# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 16:01:12 2019

@author: Richi
"""
import numpy as np

from matplotlib import cbook
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt
from matplotlib.cbook import get_sample_data
import lib.camis as camis
import csv
from scipy import ndimage

rad2deg = 180/np.pi

slopeThreshold = 30

elevationMap = np.loadtxt(open("terrainData/UMARescueAreaDEM.csv", "rb"), delimiter=" ", skiprows=0)
utmOrigin = np.loadtxt(open("terrainData/UMARescueAreaUTMOrigin.csv", "rb"), delimiter=" ", skiprows=0)
globalRes = .5


cuadrigaCamis = np.loadtxt(open("cuadriga_camis.csv", "rb"), delimiter=",", skiprows=0)

x = np.linspace(utmOrigin[0], utmOrigin[0] + globalRes*elevationMap.shape[1], elevationMap.shape[1])
y = np.linspace(utmOrigin[1], utmOrigin[1] + globalRes*elevationMap.shape[0], elevationMap.shape[0])

XX,YY = np.meshgrid(x,y)

dX,dY = np.gradient(elevationMap,*[globalRes, globalRes])

slopeMap = np.zeros_like(elevationMap)
obstacleMap = np.zeros_like(elevationMap)



for i in range(elevationMap.shape[1]):
    for j in range(elevationMap.shape[0]):
        slopeMap[j][i] = rad2deg*np.arctan(np.sqrt(dX[j][i]**2+dY[j][i]**2))
        if slopeMap[j][i] > slopeThreshold:
            obstacleMap[j][i] = 1

proximityMap = globalRes*ndimage.morphology.distance_transform_edt(1-obstacleMap)

vectorialCostMap = camis.getVectorialCostMap(slopeMap,cuadrigaCamis)

II,JJ = np.meshgrid(np.linspace(0,elevationMap.shape[1],elevationMap.shape[1]),np.linspace(0,elevationMap.shape[0],elevationMap.shape[0]))


Cad = np.sqrt(vectorialCostMap[0][:][:])

z = elevationMap
dx = globalRes
dy = globalRes


# Shade from the northwest, with the sun 45 degrees from horizontal
ls = LightSource(azdeg=315, altdeg=45)
cmap = plt.cm.gist_earth

fig, axes = plt.subplots(constrained_layout=True)

axes.contourf(II, JJ, proximityMap,100,cmap=cmap,dmax=10)
axes.set_aspect('equal')
plt.show()




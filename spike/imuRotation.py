# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 10:32:44 2020

@author: Richi
"""

from scipy.spatial.transform import Rotation as R
import numpy as np

import data.cuadrigaData.cuadriga_reader as cr
from context import camis
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yaml
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import copy
import csv
import statistics


import scipy.signal
from scipy.optimize import curve_fit

file1 = '2020_11_06_09_51_21_en_contra_de_las_agujas_del_reloj.txt'
file2 = '2020_11_06_09_49_48_a favor_de_las agujas _del_reloj.txt'
Roll = []
Pitch = []
Yaw = []

with open(file1) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=' ')
    line_count = 0
    for row in csv_reader:
        if (float(row[0]) < 0 or float(row[1]) != 0):
            Roll.append(float(row[6]))
            Pitch.append(float(row[7]))
            Yaw.append(float(row[8]))
            line_count += 1
        print(f'Processed {line_count} lines.')
#Roll = Roll[:150]
#Pitch = Pitch[:150]
#Yaw = Yaw[:150]
with open(file2) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=' ')
    line_count = 0
    for row in csv_reader:
        if float(row[0]) != 0 or float(row[1]) != 0:
            Roll.append(float(row[6]))
            Pitch.append(float(row[7]))
            Yaw.append(float(row[8]))
            line_count += 1
        print(f'Processed {line_count} lines.')
Roll = np.asarray(Roll)
Pitch = np.asarray(Pitch)
Yaw = np.asarray(Yaw)

vec = [0,0,0]
Ximu = []
Yimu = []
Zimu = []
for i,x in enumerate(Roll):
    rot = R.from_euler('ZYX', [Yaw[i], Pitch[i], Roll[i]], degrees=True)
    Ximu.append(rot.apply(np.array([1, 0, 0]),inverse = False))
    Yimu.append(rot.apply(np.array([0, 1, 0]),inverse = False))
    Zimu.append(rot.apply(np.array([0, 0, 1]),inverse = False))
#    print(vec)




#heading = np.arctan2(vy,vx)*180.0/3.1416

fig, ax = plt.subplots(figsize=(6,6), constrained_layout=True)
ax.plot(Yaw)
#ax.scatter(Yaw,[x[2] for i,x in enumerate(Ximu)])
#ax.scatter(Yaw,heading)
#ax.set_xlabel('Yaw')
#ax.set_ylabel('Slope')
#ax.plot(vx)
#ax.scatter(Yaw,Roll)
#ax.scatter(Yaw,Pitch)
#ax.scatter(Yaw,vz)

fig = plt.figure()
ax = fig.gca(projection='3d')
#ax.quiver(np.zeros(len(Ximu)), np.zeros(len(Ximu)), np.zeros(len(Ximu)), \
#          vx, vy, vz, arrow_length_ratio = 0.01,alpha=.2)
ax.quiver(0.0, 0.0, 0.0, \
          0.0, 0.0, 1.0, arrow_length_ratio = 0.05, color='c',linewidth = 5.0)
ax.quiver(0.0, 0.0, 0.0, \
          0.0, -1.0, 0.0, arrow_length_ratio = 0.05, color='b',linewidth = 5.0)
ax.scatter([x[0] for i,x in enumerate(Ximu)],\
           [x[1] for i,x in enumerate(Ximu)],\
           [x[2] for i,x in enumerate(Ximu)],c='r')
ax.scatter([x[0] for i,x in enumerate(Yimu)],\
           [x[1] for i,x in enumerate(Yimu)],\
           [x[2] for i,x in enumerate(Yimu)],c='orange')
ax.scatter([x[0] for i,x in enumerate(Zimu)],\
           [x[1] for i,x in enumerate(Zimu)],\
           [x[2] for i,x in enumerate(Zimu)],c='m')
ax.scatter(-0.5,-0.5,0.0,c='g',alpha = 0.0)
ax.scatter(0.5,0.5,0.0,c='g',alpha = 0.0)
ax.scatter(0.0,0.0,1.0,c='g')


vx = [x[0] for i,x in enumerate(Zimu)]
vy = [x[1] for i,x in enumerate(Zimu)]
vz = [x[2] for i,x in enumerate(Zimu)]
vx = np.asarray(vx)
vy = np.asarray(vy)
vz = np.asarray(vz)

vmean = [statistics.mean(vx),statistics.mean(vy),statistics.mean(vz)]
vmean = vmean / np.linalg.norm(vmean)

yawMean = np.arctan2(vmean[0],vmean[1])*180.0/3.1416
pitchMean = np.arccos(vmean[2])*180.0/3.1416

rotMean = R.from_euler('ZYX', [-yawMean, -pitchMean, 0.0], degrees=True) * R.from_euler('Z', yawMean, degrees=True)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.quiver(0.0, 0.0, 0.0, \
          0.0, 0.0, 1.0, arrow_length_ratio = 0.05, color='orange',linewidth = 5.0)
ax.quiver(0.0, 0.0, 0.0, \
          vmean[0],vmean[1],vmean[2], arrow_length_ratio = 0.05, color='b',linewidth = 5.0)
ax.quiver(0.0, 0.0, 0.0, \
          0.0, -1.0, 0.0, arrow_length_ratio = 0.05, color='m',linewidth = 5.0)
#ax.scatter(vx,vy,vz,c='r')
ax.scatter(-0.5,-0.5,0.0,c='g',alpha = 0.0)
ax.scatter(0.5,0.5,0.0,c='g',alpha = 0.0)
ax.scatter(0.0,0.0,1.0,c='g')

#yawmean = np.arctan2(vmean[1],vmean[0])*180.0/3.1416
#pitchmean = - np.arccos(vmean[2])*180.0/3.1416

pitchdef = - np.arctan2(.1,np.sqrt(1-.1**2))*180.0/3.1416

#zvector = np.asarray([0.0,0.0,1.0])
#ResultingRotMatrix = R.align_vectors(zvector.transpose(),vmean.transpose())

#index 549 -> 0 degrees
#index 102 -> 80 degrees
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.quiver(0.0, 0.0, 0.0, \
          1.0, 0.0, 0.0, arrow_length_ratio = 0.05, color='orange',linewidth = 5.0)
ax.quiver(0.0, 0.0, 0.0, \
          0.0, 1.0, 0.0, arrow_length_ratio = 0.05, color='orange',linewidth = 5.0)
ax.quiver(0.0, 0.0, 0.0, \
          0.0, 0.0, 1.0, arrow_length_ratio = 0.05, color='orange',linewidth = 5.0)
ax.quiver(0.0, 0.0, 0.0, \
          Ximu[102][0], Ximu[102][1], Ximu[102][2], arrow_length_ratio = 0.05, color='m',linewidth = 5.0)
ax.quiver(0.0, 0.0, 0.0, \
          Yimu[102][0], Yimu[102][1], Yimu[102][2], arrow_length_ratio = 0.05, color='m',linewidth = 5.0)
ax.quiver(0.0, 0.0, 0.0, \
          Zimu[102][0], Zimu[102][1], Zimu[102][2], arrow_length_ratio = 0.05, color='m',linewidth = 5.0)
ax.scatter(-1.0,-1.0,0.0,c='g',alpha = 0.0)
ax.scatter(1.0,1.0,1.0,c='g',alpha = 0.0)
ax.scatter(0.0,0.0,1.0,c='g')

Rimu = (
[
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
])
    
Rcuad = (
[
    Ximu[549],
    Yimu[549],
    Zimu[549],
])

I2Cmatrix = R.from_matrix(Rcuad)
    
#I2Cmatrix,_ = R.align_vectors(Rcuad,Rimu)
I2Cmatrix.as_euler('ZYX',degrees=True)
#I2Cmatrix = R.from_euler('Z',45.0,degrees=True)*I2Cmatrix
#I2Cmatrix = I2Cmatrix.inv()
I2Cmatrix.as_euler('ZYX',degrees=True)

I2Cmatrix = I2Cmatrix.as_matrix()

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.quiver(0.0, 0.0, 0.0, \
          1.0, 0.0, 0.0, arrow_length_ratio = 0.05, color='orange',linewidth = 5.0)
ax.quiver(0.0, 0.0, 0.0, \
          0.0, 1.0, 0.0, arrow_length_ratio = 0.05, color='orange',linewidth = 5.0)
ax.quiver(0.0, 0.0, 0.0, \
          0.0, 0.0, 1.0, arrow_length_ratio = 0.05, color='orange',linewidth = 5.0)
ax.quiver(0.0, 0.0, 0.0, \
          I2Cmatrix[0][0], I2Cmatrix[0][1], I2Cmatrix[0][2], arrow_length_ratio = 0.05, color='m',linewidth = 5.0)
ax.quiver(0.0, 0.0, 0.0, \
          I2Cmatrix[1][0], I2Cmatrix[1][1], I2Cmatrix[1][2], arrow_length_ratio = 0.05, color='m',linewidth = 5.0)
ax.quiver(0.0, 0.0, 0.0, \
          I2Cmatrix[2][0], I2Cmatrix[2][1], I2Cmatrix[2][2], arrow_length_ratio = 0.05, color='m',linewidth = 5.0)
ax.scatter(-1.0,-1.0,0.0,c='g',alpha = 0.0)
ax.scatter(1.0,1.0,1.0,c='g',alpha = 0.0)
ax.scatter(0.0,0.0,1.0,c='g')


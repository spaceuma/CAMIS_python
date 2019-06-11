# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 12:44:19 2019

@author: Richi
"""
import csv
import utm
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal

deg2rad = np.pi/180
rad2deg = 180/np.pi

def readCuadrigaData(file):
    Roll = []
    Pitch = []
    Yaw = []
    Current = []
    utmPoseX = []
    utmPoseY = []
    Time = []
    
    with open(file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        line_count = 0
        for row in csv_reader:
            error = float(row[14])
            reference = float(row[0])
            if error < 0.3 and reference > 0:
                c = np.abs(float(row[4])) + np.abs(float(row[5]))
                r = -float(row[6])-5.5
                p = float(row[7])+13.7
                utmData = utm.from_latlon(float(row[9]), float(row[10]))
                utmPoseX.append(utmData[0])
                utmPoseY.append(utmData[1])
                Roll.append(r)
                Pitch.append(p)
                Current.append(c)
                t = dt.datetime.strptime(row[15], '%H:%M:%S.%f')
                Time.append((t.hour * 60 + t.minute) * 60 + t.second + t.microsecond/1000000)
            line_count += 1
        print(f'Processed {line_count} lines.')
        
    # Adjustments
    initTime = Time[0]
    Time[:] = [t - initTime for t in Time]
    
    #Smooth the path
    utmPoseX = scipy.ndimage.median_filter(utmPoseX,size=5)
    utmPoseY = scipy.ndimage.median_filter(utmPoseY,size=5)
    
    #Compute Yaw
    dX = np.append(0,np.diff(utmPoseX))
    dY = np.append(0,np.diff(utmPoseY))
    Yaw = []
    for i,y in enumerate(dX):
        Yaw.append(rad2deg*np.arctan2(dY[i],dX[i]))
    
    #Eliminate unwanted values
    Roll[:] = [r for i,r in enumerate(Roll) if (dX[i] != 0 or dY[i] != 0) and (Current[i] > 4) and (Current[i]) < 40]
    Pitch[:] = [p for i,p in enumerate(Pitch) if (dX[i] != 0 or dY[i] != 0) and (Current[i] > 4) and (Current[i]) < 40]
    Yaw[:] = [y for i,y in enumerate(Yaw) if (dX[i] != 0 or dY[i] != 0) and (Current[i] > 4) and (Current[i]) < 40]
    Current[:] = [c for i,c in enumerate(Current) if (dX[i] != 0 or dY[i] != 0) and (Current[i] > 4) and (Current[i]) < 40]
    
    # Filtering
    Current = scipy.ndimage.median_filter(Current,size=10)
    Roll = scipy.ndimage.median_filter(Roll,size=10)
    Pitch = scipy.ndimage.median_filter(Pitch,size=10)
            
    return Roll, Pitch, Yaw, np.ndarray.tolist(Current)


def getData(file):
    Roll = []
    Pitch = []
    Yaw = []
    Slope = []
    Current = []
    utmPoseX = []
    utmPoseY = []
    utmPose = []
    GPSHeading = []
    
    Time = []
    
    deg2rad = np.pi/180
    rad2deg = 180/np.pi
    
    with open(file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        line_count = 0
        for row in csv_reader:
            error = float(row[14])
            reference = float(row[0])
            if error < 0.3 and reference > 0:
                c = np.abs(float(row[4])) + np.abs(float(row[5]))
                r = -float(row[6])-5.5
                p = float(row[7])+13.7
                utmData = utm.from_latlon(float(row[9]), float(row[10]))
                utmPoseX.append(utmData[0])
                utmPoseY.append(utmData[1])
                utmPose.append([utmData[0], utmData[1]])
                Roll.append(r)
                Pitch.append(p)
                Yaw.append(-float(row[8])-50.0)
                Current.append(c)
                GPSHeading.append(-float(row[12])+90.0)
                t = dt.datetime.strptime(row[15], '%H:%M:%S.%f')
                Time.append((t.hour * 60 + t.minute) * 60 + t.second + t.microsecond/1000000)
            line_count += 1
        print(f'Processed {line_count} lines.')
        
    # Adjustments
    initTime = Time[0]
    Time[:] = [t - initTime for t in Time]
    
    # Filtering
    Current = scipy.ndimage.median_filter(Current,size=10)
    Roll = scipy.ndimage.median_filter(Roll,size=10)
    Pitch = scipy.ndimage.median_filter(Pitch,size=10)
    GPSHeading = scipy.ndimage.median_filter(GPSHeading,size=10)
    utmPoseX = scipy.ndimage.median_filter(utmPoseX,size=5)
    utmPoseY = scipy.ndimage.median_filter(utmPoseY,size=5)
    
    # Gradient is computed
    Gradient = np.zeros_like(Roll)
    for i, r in enumerate(Roll):
        Gradient[i] = rad2deg*np.arccos(np.cos(deg2rad*Roll[i])*np.cos(deg2rad*Pitch[i]))
        
    # beta angle is computed
    BetaX = np.zeros_like(Roll)
    BetaY = np.zeros_like(Roll)
    aspectVectorX = []
    aspectVectorY = []
    headingVectorX = []
    headingVectorY=[]
    for i, r in enumerate(Roll):
        aspect = [np.cos(deg2rad*GPSHeading[i])*np.sin(deg2rad*Pitch[i])*np.cos(deg2rad*Roll[i]) + np.sin(deg2rad*GPSHeading[i])*np.sin(deg2rad*Roll[i]), np.sin(deg2rad*GPSHeading[i])*np.sin(deg2rad*Pitch[i])*np.cos(deg2rad*Roll[i]) - np.cos(deg2rad*GPSHeading[i])*np.sin(deg2rad*Roll[i])]
        aspect = aspect/np.linalg.norm(aspect)
        aspectVectorX.append(aspect[0])
        aspectVectorY.append(aspect[1])
        heading = [np.cos(deg2rad*GPSHeading[i])*np.cos(deg2rad*Pitch[i]), np.sin(deg2rad*GPSHeading[i])*np.cos(deg2rad*Pitch[i])]
        heading = heading/np.linalg.norm(heading)
        headingVectorX.append(heading[0])
        headingVectorY.append(heading[1])
        c, s = aspect[0], aspect[1]
        R = np.array(((c,s), (-s, c)))
        T = np.dot(R,np.transpose(heading))
        BetaX[i] = T[0]
        BetaY[i] = T[1]
        
    # Heading is computed
    dX = np.append(0,np.diff(utmPoseX))
    dY = np.append(0,np.diff(utmPoseY))
    Segments = []
    movingTime = []
    movingHeading = []
    movingCurrent = []
    movingGradient = []
    movingBetaX = []
    movingBetaY = []
    for i, d in enumerate(dX):
        if (dX[i] != 0 or dY[i] != 0) and (Current[i] > 4) and (Current[i]) < 40:
            Segments.append(np.sqrt(np.power(dX[i],2) + np.power(dY[i],2)))
            movingTime.append(Time[i])
            movingHeading.append(rad2deg*np.arctan2(dY[i],dX[i]))
            movingCurrent.append(Current[i])
            movingGradient.append(Gradient[i])
            movingBetaX.append(BetaX[i])
            movingBetaY.append(BetaY[i])
    return movingBetaX, movingBetaY, movingCurrent, movingGradient

def showData(file):
    

    Roll = []
    Pitch = []
    Yaw = []
    Slope = []
    Current = []
    utmPoseX = []
    utmPoseY = []
    utmPose = []
    GPSHeading = []
    
    Time = []
    
    deg2rad = np.pi/180
    rad2deg = 180/np.pi
    
    with open(file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        line_count = 0
        for row in csv_reader:
            error = float(row[14])
            reference = float(row[0])
            if error < 0.3 and reference > 0:
                c = np.abs(float(row[4])) + np.abs(float(row[5]))
                r = -float(row[6])-5.5
                p = float(row[7])+13.7
                utmData = utm.from_latlon(float(row[9]), float(row[10]))
                utmPoseX.append(utmData[0])
                utmPoseY.append(utmData[1])
                utmPose.append([utmData[0], utmData[1]])
                Roll.append(r)
                Pitch.append(p)
                Yaw.append(-float(row[8])-50.0)
                Current.append(c)
                GPSHeading.append(-float(row[12])+90.0)
                t = dt.datetime.strptime(row[15], '%H:%M:%S.%f')
                Time.append((t.hour * 60 + t.minute) * 60 + t.second + t.microsecond/1000000)
            line_count += 1
        print(f'Processed {line_count} lines.')
        
    # Adjustments
    initTime = Time[0]
    Time[:] = [t - initTime for t in Time]
    
    # Filtering
    Current = scipy.ndimage.median_filter(Current,size=10)
    Roll = scipy.ndimage.median_filter(Roll,size=10)
    Pitch = scipy.ndimage.median_filter(Pitch,size=10)
    GPSHeading = scipy.ndimage.median_filter(GPSHeading,size=10)
    utmPoseX = scipy.ndimage.median_filter(utmPoseX,size=5)
    utmPoseY = scipy.ndimage.median_filter(utmPoseY,size=5)
    
    # Gradient is computed
    Gradient = np.zeros_like(Roll)
    for i, r in enumerate(Roll):
        Gradient[i] = rad2deg*np.arccos(np.cos(deg2rad*Roll[i])*np.cos(deg2rad*Pitch[i]))
        
    # beta angle is computed
    BetaX = np.zeros_like(Roll)
    BetaY = np.zeros_like(Roll)
    aspectVectorX = []
    aspectVectorY = []
    headingVectorX = []
    headingVectorY=[]
    for i, r in enumerate(Roll):
        aspect = [np.cos(deg2rad*GPSHeading[i])*np.sin(deg2rad*Pitch[i])*np.cos(deg2rad*Roll[i]) + np.sin(deg2rad*GPSHeading[i])*np.sin(deg2rad*Roll[i]), np.sin(deg2rad*GPSHeading[i])*np.sin(deg2rad*Pitch[i])*np.cos(deg2rad*Roll[i]) - np.cos(deg2rad*GPSHeading[i])*np.sin(deg2rad*Roll[i])]
        aspect = aspect/np.linalg.norm(aspect)
        aspectVectorX.append(aspect[0])
        aspectVectorY.append(aspect[1])
        heading = [np.cos(deg2rad*GPSHeading[i])*np.cos(deg2rad*Pitch[i]), np.sin(deg2rad*GPSHeading[i])*np.cos(deg2rad*Pitch[i])]
        heading = heading/np.linalg.norm(heading)
        headingVectorX.append(heading[0])
        headingVectorY.append(heading[1])
        c, s = aspect[0], aspect[1]
        R = np.array(((c,s), (-s, c)))
        T = np.dot(R,np.transpose(heading))
        BetaX[i] = T[0]
        BetaY[i] = T[1]
        
    # Heading is computed
    dX = np.append(0,np.diff(utmPoseX))
    dY = np.append(0,np.diff(utmPoseY))
    Segments = []
    movingTime = []
    movingHeading = []
    movingCurrent = []
    movingGradient = []
    movingBetaX = []
    movingBetaY = []
    for i, d in enumerate(dX):
        if (dX[i] != 0 or dY[i] != 0) and (Current[i] > 4) and (Current[i]) < 40:
            Segments.append(np.sqrt(np.power(dX[i],2) + np.power(dY[i],2)))
            movingTime.append(Time[i])
            movingHeading.append(rad2deg*np.arctan2(dY[i],dX[i]))
            movingCurrent.append(Current[i])
            movingGradient.append(Gradient[i])
            movingBetaX.append(BetaX[i])
            movingBetaY.append(BetaY[i])
    
    # Energy Consumed array
    Energy = np.zeros_like(Current)
    for i, element in enumerate(Energy):
        Energy[i] = np.trapz(Current[0:i+1],Time[0:i+1])
        
    fig, ax = plt.subplots()
    ax.plot(Time,Roll,Time,Pitch,Time,Gradient)
    ax.set_xlim(0, Time[-1])
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Angle [degrees]')
    ax.grid(True)
    ax.legend(labels = ['Roll','Pitch','Gradient'])
    
    fig, ax = plt.subplots()
    ax.plot(Time,GPSHeading,Time,Yaw,movingTime,movingHeading)
    ax.set_xlim(0, Time[-1])
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Orientation [degrees]')
    ax.grid(True)
    ax.legend(labels = ['GPS Heading','Yaw','Computed Heading'])
    
    
    fig, ax = plt.subplots()
    ax.plot(movingTime,movingCurrent)
    ax.set_xlim(0, Time[-1])
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Current [A]')
    ax.grid(True)
    ax.legend(labels = ['Current','filtered'])
    
    fig, ax = plt.subplots()
    ax.plot(utmPoseX,utmPoseY)
    ax.set_xlabel('UTM-X [m]')
    ax.set_ylabel('UTM-Y [m]')
    ax.grid(True)
    ax.legend(labels = ['Path'])
    plt.tight_layout()
    ax.set_aspect('equal')
    Q1 = ax.quiver(utmPoseX, utmPoseY, headingVectorX, headingVectorY,units='width',scale=50,color = 'g')
    ax.quiver(utmPoseX, utmPoseY, aspectVectorX, aspectVectorY,units='width',scale=50,color = 'r')
    ax.quiverkey(Q1, X=0.3, Y=1.0, U=10,
                 label='Quiver key, length = 10', labelpos='E')
    
    fig, ax = plt.subplots()
    ax.stackplot(Time,Energy)
    ax.set_xlim(0, Time[-1])
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Energy [As]')
    ax.grid(True)
    ax.legend(labels = ['Energy consumed'])
    
#showData('20190531/JornadasRescate01.txt')
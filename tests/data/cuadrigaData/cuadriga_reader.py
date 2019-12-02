# -*- coding: utf-8 -*-
"""
Library to extract data from Cuadriga
"""
import csv
import utm
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal

deg2rad = np.pi/180
rad2deg = 180/np.pi

def readCuadrigaData(file, show=0):
    Roll = []
    Pitch = []
    Yaw = []
    Current = []
    Speed = []
    utmPoseX = []
    utmPoseY = []
    Time = []
    Ref = []
    Error = []
    GPSspeed = []
    
    with open(file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        line_count = 0
        for row in csv_reader:
            error = float(row[14])
            reference = float(row[0])
            if error < 1.0 and reference > 0:
                c = np.abs(float(row[4])) + np.abs(float(row[5]))
                r = -float(row[6])-5.5
#                p = float(row[7])+13.7
                p = float(row[7])+10.7
                utmData = utm.from_latlon(float(row[9]), float(row[10]))
                utmPoseX.append(utmData[0])
                utmPoseY.append(utmData[1])
                GPSspeed.append(float(row[13]))
                Roll.append(r)
                Pitch.append(p)
                Current.append(c)
                Ref.append(reference)
                Error.append(error)
                t = dt.datetime.strptime(row[15], '%H:%M:%S.%f')
                Time.append((t.hour * 60 + t.minute) * 60 + t.second + t.microsecond/1000000)
            line_count += 1
        print(f'Processed {line_count} lines.')
    
    #Smooth the path
    N = 10
    utmPoseX = np.convolve(utmPoseX, np.ones((N,))/N, mode='same')
    utmPoseX = utmPoseX[5:-5]
    utmPoseY = np.convolve(utmPoseY, np.ones((N,))/N, mode='same')
    utmPoseY = utmPoseY[5:-5]
#    Current = scipy.ndimage.median_filter(Current,size=10)
    Current = np.convolve(Current, np.ones((N,))/N, mode='same')
    Roll = np.convolve(Roll, np.ones((N,))/N, mode='same')
    Pitch = np.convolve(Pitch, np.ones((N,))/N, mode='same')
    Time = Time[5:-5]
    Roll = Roll[5:-5]
    Pitch = Pitch[5:-5]
    Current = Current[5:-5]
    GPSspeed = GPSspeed[5:-5]
    
    # Adjustments
    initTime = Time[0]
    Time[:] = [t - initTime for t in Time]
    dX = np.append(0,np.diff(utmPoseX))
    dY = np.append(0,np.diff(utmPoseY))
    segmentPath = np.sqrt(dX**2+dY**2)
    dT = np.append(0,np.diff(Time))
    Speed = segmentPath/dT
    
    logic = np.logical_and((segmentPath > 0),(dT < 0.15))
    logic = np.logical_and((segmentPath < 0.075), logic)
    logic = np.logical_and((dT > 0.0), logic)
    
    Time[:] = [t for i,t in enumerate(Time) if logic[i]]
    utmPoseX = [x for i,x in enumerate(utmPoseX) if logic[i]]
    utmPoseY = [y for i,y in enumerate(utmPoseY) if logic[i]]
    Current = [y for i,y in enumerate(Current) if logic[i]]
    dX = [x for i,x in enumerate(dX) if logic[i]]
    dY = [y for i,y in enumerate(dY) if logic[i]]
    dT = [y for i,y in enumerate(dT) if logic[i]]
    segmentPath = [y for i,y in enumerate(segmentPath) if logic[i]]
    Roll = [y for i,y in enumerate(Roll) if logic[i]]
    Pitch = [y for i,y in enumerate(Pitch) if logic[i]]
    Yaw = [y for i,y in enumerate(Yaw) if logic[i]]
    GPSspeed = [s for i,s in enumerate(GPSspeed) if logic[i]]
    

    heading = np.arctan2(dY,dX)
    dHeading = np.zeros_like(heading)
    for i in range(1,len(heading)):
        dHeading[i] = np.arccos((dX[i]*dX[i-1]+dY[i]*dY[i-1])/(segmentPath[i]*segmentPath[i-1]))
        if np.isnan(dHeading[i]):
            dHeading[i] = 0.0
    
    angularSpeed = np.asarray(dHeading)/np.asarray(dT)

    Speed = np.asarray(segmentPath)/np.asarray(dT)
    
    Yaw = rad2deg*np.arctan2(dY,dX)
    
    Curvature = np.abs(angularSpeed/Speed)
    Curvature = np.convolve(Curvature, np.ones((N,))/N, mode='same')
    
    logic = (Curvature < 0.45)
    
    utmPoseX = [x for i,x in enumerate(utmPoseX) if logic[i]]
    utmPoseY = [y for i,y in enumerate(utmPoseY) if logic[i]]
    Current = [y for i,y in enumerate(Current) if logic[i]]
    Time = [t for i,t in enumerate(Time) if logic[i]]
    heading = [y for i,y in enumerate(heading) if logic[i]]
    dHeading = [y for i,y in enumerate(dHeading) if logic[i]]
    segmentPath = [y for i,y in enumerate(segmentPath) if logic[i]]
    dT = [y for i,y in enumerate(dT) if logic[i]]
    Roll = [y for i,y in enumerate(Roll) if logic[i]]
    Pitch = [y for i,y in enumerate(Pitch) if logic[i]]
    Yaw = [y for i,y in enumerate(Yaw) if logic[i]]
    GPSspeed = [s for i,s in enumerate(GPSspeed) if logic[i]]
    Speed = [s for i,s in enumerate(Speed) if logic[i]]
    Curvature = [y for i,y in enumerate(Curvature) if logic[i]]
    
    
    Time = Time[5:]
    dHeading = dHeading[5:]
    segmentPath = segmentPath[5:]
    traversedDist = np.cumsum(segmentPath)
    dT = dT[5:]
    Curvature = Curvature[5:]
    
    
    
    if show == 1:
        fig, ax = plt.subplots()
        ax.plot(Time,Speed,Time,GPSspeed)
#        ax.set_xlim(0, Time[-1])
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Speed [m/s]')
        ax.grid(True)
        ax.legend(labels = ['Computed Speed','GPS Speed'])
        
        fig, ax = plt.subplots()
        ax.plot(utmPoseX,utmPoseY)
        ax.set_xlabel('UTM-X [m]')
        ax.set_ylabel('UTM-Y [m]')
        ax.grid(True)
        ax.legend(labels = ['Path'])
        plt.tight_layout()
        ax.set_aspect('equal')
        
        Gradient = np.zeros_like(Roll)
        for i, r in enumerate(Roll):
            Gradient[i] = rad2deg*np.arccos(np.cos(deg2rad*Roll[i])*np.cos(deg2rad*Pitch[i]))
        
        fig, ax = plt.subplots()
        ax.plot(Time,Roll,Time,Pitch,Time,Gradient)
        ax.set_xlim(0, Time[-1])
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Angle [degrees]')
        ax.grid(True)
        ax.legend(labels = ['Roll','Pitch','Gradient'])
        
        fig, ax = plt.subplots()
        ax.plot(Time,Current)
#        ax.set_xlim(0, Time[-1])
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('[A]')
        ax.grid(True)
        ax.legend(labels = ['Current'])
        
        fig, ax = plt.subplots()
        ax.plot(Time,Curvature)
#        ax.set_xlim(0, Time[-1])
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('[A]')
        ax.grid(True)
        ax.legend(labels = ['Curvature'])
            
    return utmPoseX[5:], utmPoseY[5:], heading[5:], Roll[5:], Pitch[5:], Yaw[5:], Current[5:], Speed[5:], traversedDist


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
    ax.plot(Time,Current,movingTime,movingCurrent)
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
    
def getCuadrigaPath(file):
    Roll = []
    Pitch = []
    Yaw = []
    Slope = []
    Current = []
    utmPoseX = []
    utmPoseY = []
    utmPose = []
    GPSHeading = []
    dHeading = []
    
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
    #Smooth the path
    N = 10
    utmPoseX = np.convolve(utmPoseX, np.ones((N,))/N, mode='same')
    utmPoseX = utmPoseX[5:-5]
    utmPoseY = np.convolve(utmPoseY, np.ones((N,))/N, mode='same')
    utmPoseY = utmPoseY[5:-5]
    Current = scipy.ndimage.median_filter(Current,size=10)
    Time = Time[5:-5]
    Roll = Roll[5:-5]
    Pitch = Pitch[5:-5]
    Current = Current[5:-5]
        
    # Adjustments
    initTime = Time[0]
    Time[:] = [t - initTime for t in Time]
    dX = np.append(0,np.diff(utmPoseX))
    dY = np.append(0,np.diff(utmPoseY))
    segmentPath = np.sqrt(dX**2+dY**2)
    dT = np.append(0,np.diff(Time))
    Speed = segmentPath/dT
    
    logic = np.logical_and((segmentPath > 0),(dT < 0.15))
    logic = np.logical_and((segmentPath < 0.075), logic)
    logic = np.logical_and((dT > 0.0), logic)
    
    Time[:] = [t for i,t in enumerate(Time) if logic[i]]
    utmPoseX = [x for i,x in enumerate(utmPoseX) if logic[i]]
    utmPoseY = [y for i,y in enumerate(utmPoseY) if logic[i]]
    Current = [y for i,y in enumerate(Current) if logic[i]]
    dX = [x for i,x in enumerate(dX) if logic[i]]
    dY = [y for i,y in enumerate(dY) if logic[i]]
    dT = [y for i,y in enumerate(dT) if logic[i]]
    segmentPath = [y for i,y in enumerate(segmentPath) if logic[i]]

    
    
    heading = np.arctan2(dY,dX)
    dHeading[:] = np.zeros_like(heading)
    for i in range(1,len(heading)):
        dHeading[i] = np.arccos((dX[i]*dX[i-1]+dY[i]*dY[i-1])/(segmentPath[i]*segmentPath[i-1]))
        if np.isnan(dHeading[i]):
            dHeading[i] = 0.0
    
    angularSpeed = np.asarray(dHeading)/np.asarray(dT)
    
    
    
    Speed = np.asarray(segmentPath)/np.asarray(dT)
    
    Curvature = np.abs(angularSpeed/Speed)
    Curvature = np.convolve(Curvature, np.ones((N,))/N, mode='same')
    
    logic = (Curvature < 0.45)
    
    utmPoseX = [x for i,x in enumerate(utmPoseX) if logic[i]]
    utmPoseY = [y for i,y in enumerate(utmPoseY) if logic[i]]
    Current = [y for i,y in enumerate(Current) if logic[i]]
    Time = [t for i,t in enumerate(Time) if logic[i]]
    dHeading = [y for i,y in enumerate(dHeading) if logic[i]]
    segmentPath = [y for i,y in enumerate(segmentPath) if logic[i]]
    dT = [y for i,y in enumerate(dT) if logic[i]]
    Curvature = [y for i,y in enumerate(Curvature) if logic[i]]
    
    utmPoseX = utmPoseX[5:]
    utmPoseY = utmPoseY[5:]
    Current = Current[5:]
    Time = Time[5:]
    dHeading = dHeading[5:]
    segmentPath = segmentPath[5:]
    dT = dT[5:]
    Curvature = Curvature[5:]
    
    return np.asarray(utmPoseX), np.asarray(utmPoseY), np.asarray(Current),\
           np.asarray(Curvature), np.asarray(dHeading),\
           np.asarray(segmentPath), np.asarray(dT)
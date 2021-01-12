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
from scipy.spatial.transform import Rotation as R

deg2rad = np.pi/180
rad2deg = 180/np.pi

def readCuadrigaData(file, counter_lim):
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
    Slope = []
    
    reference = 0
    
    path_counter = 0
    
    I2C = R.from_euler('ZYX', [159.68348098,  14.18826271,   5.65423381], degrees=True)
#    I2C = R.from_euler('ZYX', [0.0, -11.4, -10.2 ], degrees=True)
    I2Cinv = I2C.inv()
#    horizontalRot = R.from_euler('zyx', [0.0, -11.4, -10.1 ], degrees=True)
#    horizontalRot = R.from_euler('zyx', [-30.0, -13.7, -5.5 ], degrees=True)
#    horizontalRot = R.from_euler('zyx', [-60.3, -15.3, 1.1 ], degrees=True)
#    horizontalRot = R.from_euler('zyx', [-90.0, -8.2, 6.2 ], degrees=True)
#    hR = R.from_euler('z', 0.0, degrees=True) * horizontalRot.inv()  
    
#    I2C_Z = [-0.1377589 ,  0.15973515,  0.97750047]
#    I2C_Z = [-0.19453349,  0.17708474,  0.96477858]
    I2C_Z = [0.0,0.0,0.0]
    with open(file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        line_count = 0
        for row in csv_reader:
            error = float(row[14])
            if float(row[0]) > 0 and reference == 0:
                path_counter = path_counter + 1
            reference = float(row[0])
            if error < 1.0 and reference > 0 and path_counter == counter_lim:
                imuRot = R.from_euler('ZYX', [float(row[8]), float(row[7]), float(row[6])], degrees=True)
                #Z component of the orientation vector from IMU perspective
                imuRot = imuRot.as_matrix()
                imuZ = [imuRot[0][2], imuRot[1][2], imuRot[2][2]]
                if np.linalg.norm(I2C_Z) < 0.1:
                    I2C_Z = imuZ
                Slope.append(np.arccos(np.sum(np.dot(I2C_Z,imuZ)))*180.0/3.1416)
                c = np.abs(float(row[4])) + np.abs(float(row[5]))
                I2R = R.from_euler('ZYX', [float(row[8]), float(row[7]), float(row[6])], degrees=True)
#                fixRot = R.from_euler('XYZ', [-float(row[6]), 12.177, - float(row[8])], degrees=True)
                C2R = I2Cinv*I2R
                R2C = C2R.inv()
                y,p,r = R2C.as_euler('ZYX', degrees=True)
#                r = float(row[6]) + 6.0
##                r = -float(row[6])-5.5
##                p = float(row[7])+13.7
##                p = float(row[7])+10.7
#                p = float(row[7])+17
#                y = float(row[8])
                utmData = utm.from_latlon(float(row[9]), float(row[10]))
                utmPoseX.append(utmData[0])
                utmPoseY.append(utmData[1])
                GPSspeed.append(float(row[13]))
                Roll.append(float(row[6]))
                Pitch.append(float(row[7]))
                Yaw.append(float(row[8]))
                Current.append(c)
                Ref.append(reference)
                Error.append(error)
                t = dt.datetime.strptime(row[15], '%H:%M:%S.%f')
                Time.append((t.hour * 60 + t.minute) * 60 + t.second + t.microsecond/1000000)
            if reference == 0 and path_counter == counter_lim:
                break
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
    Yaw = Yaw[5:-5]
    Slope = Slope[5:-5]
    Current = Current[5:-5]
    GPSspeed = GPSspeed[5:-5]
    
    # Adjustments
    initTime = Time[0]
    Time[:] = [ti - initTime for ti in Time]
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
    Slope = [y for i,y in enumerate(Slope) if logic[i]]
    GPSspeed = [s for i,s in enumerate(GPSspeed) if logic[i]]
    

    heading = np.arctan2(dY,dX)
    dHeading = np.zeros_like(heading)
    for i in range(1,len(heading)):
        dHeading[i] = np.arccos((dX[i]*dX[i-1]+dY[i]*dY[i-1])/(segmentPath[i]*segmentPath[i-1]))
        if np.isnan(dHeading[i]):
            dHeading[i] = 0.0
    
    angularSpeed = np.asarray(dHeading)/np.asarray(dT)

    Speed = np.asarray(segmentPath)/np.asarray(dT)
    
#    Yaw = rad2deg*np.arctan2(dY,dX)
    
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
    Slope = [y for i,y in enumerate(Slope) if logic[i]]
    GPSspeed = [s for i,s in enumerate(GPSspeed) if logic[i]]
    Speed = [s for i,s in enumerate(Speed) if logic[i]]
    Curvature = [y for i,y in enumerate(Curvature) if logic[i]]
    
    
    Time = Time[5:]
    Time = [t - Time[0] for i,t in enumerate(Time)]
    dHeading = dHeading[5:]
    segmentPath = segmentPath[5:]
    traversedDist = np.cumsum(segmentPath)
    dT = dT[5:]
    Curvature = Curvature[5:]
    Speed = Speed[5:]
    Current = Current[5:]
    GPSspeed = GPSspeed[5:]
    Roll = Roll[5:]
    Pitch = Pitch[5:]
    Yaw = Yaw[5:]
    utmPoseX = utmPoseX[5:]
    utmPoseY = utmPoseY[5:]
    heading = heading[5:]
    Slope = Slope[5:]
    
    segmentPath = np.zeros_like(utmPoseX)
    for index, waypoint in enumerate(utmPoseX):
            if index == 0:
                segmentPath[index] = 0.0
            else:
                A = [utmPoseX[index] - utmPoseX[index-1],utmPoseY[index] - utmPoseY[index-1]]
                segmentPath[index] = np.linalg.norm(A)               
    traversedDist = np.cumsum(segmentPath)
    
    
#    Rimu = (
#    [
#        [1, 0, 0],
#        [0, 1, 0],
#        [0, 0, 1],
#    ])
#        
#    Rcuad = (
#    [
#        [0.98027117, 0.        , 0.19765734],
#        [ 0.0350021 ,  0.98419561, -0.17359107],
#        [-0.19453349,  0.17708474,  0.96477858],
#    ])
    
##    I2Cmatrix,_ = R.align_vectors(Rcuad,Rimu)
##    I2Cmatrix = R.from_matrix(Rcuad)
#    # Make here a more sofisticated calibration
#    I2Cmatrix = R.from_euler('ZYX', [Yaw[0], Pitch[0], Roll[0]], degrees=True)
#    I2Cmatrix.as_euler('ZYX',degrees=True)
##    Y0 = Yaw[0]
#    for i,x in enumerate(Roll):
#        if i == 0:
#            A = [utmPoseX[1] - utmPoseX[0], utmPoseY[1] - utmPoseY[0]]
#            heading[i] = np.arctan2(A[1],A[0]) * 180.0/3.1416
#        elif i == len(utmPoseX) - 1:
#            A = [utmPoseX[i] - utmPoseX[i-1], utmPoseY[i] - utmPoseY[i-1]]
#            heading[i] = np.arctan2(A[1],A[0]) * 180.0/3.1416
#        else:
#            A1 = [utmPoseX[i] - utmPoseX[i-1], utmPoseY[i] - utmPoseY[i-1]]
#            A2 = [utmPoseX[i+1] - utmPoseX[i], utmPoseY[i+1] - utmPoseY[i]]
#            heading[i] = np.arctan2(A1[1]+A2[1],A1[0]+A1[0]) * 180.0/3.1416
##        I2C = I2Cmatrix*R.from_euler('z',Yaw[i] + heading[i] + 0.0,degrees=True)
##        I2C = I2Cmatrix*R.from_euler('z',-Yaw[i]+Y0,degrees=True)
#        I2C = I2Cmatrix
#        I2Cinv = I2C.inv()
#        I2R = R.from_euler('ZYX', [Yaw[i], Pitch[i], Roll[i]], degrees=True)
#        C2R = I2C.inv()*I2R
##        R2C = C2R.inv()
#        Yaw[i],Pitch[i],Roll[i] = C2R.as_euler('ZYX', degrees=True)
    
#    if show == 1:
#        fig, ax = plt.subplots()
#        ax.plot(Time,Speed,Time,GPSspeed)
##        ax.set_xlim(0, Time[-1])
#        ax.set_xlabel('Time [s]')
#        ax.set_ylabel('Speed [m/s]')
#        ax.grid(True)
#        ax.legend(labels = ['Computed Speed','GPS Speed'])
#        
#        fig, ax = plt.subplots()
#        ax.plot(utmPoseX,utmPoseY)
#        ax.set_xlabel('UTM-X [m]')
#        ax.set_ylabel('UTM-Y [m]')
#        ax.grid(True)
#        ax.legend(labels = ['Path'])
#        plt.tight_layout()
#        ax.set_aspect('equal')
#        
#        Gradient = np.zeros_like(Roll)
#        for i, r in enumerate(Roll):
#            Gradient[i] = rad2deg*np.arccos(np.cos(deg2rad*Roll[i])*np.cos(deg2rad*Pitch[i]))
#        
#        fig, ax = plt.subplots()
#        ax.plot(Time,Roll,Time,Pitch,Time,Gradient)
#        ax.set_xlim(0, Time[-1])
#        ax.set_xlabel('Time [s]')
#        ax.set_ylabel('Angle [degrees]')
#        ax.grid(True)
#        ax.legend(labels = ['Roll','Pitch','Gradient'])
#        
#        fig, ax = plt.subplots()
#        ax.plot(Time,Current)
##        ax.set_xlim(0, Time[-1])
#        ax.set_xlabel('Time [s]')
#        ax.set_ylabel('[A]')
#        ax.grid(True)
#        ax.legend(labels = ['Current'])
#        
#        fig, ax = plt.subplots()
#        ax.plot(Time,Curvature)
##        ax.set_xlim(0, Time[-1])
#        ax.set_xlabel('Time [s]')
#        ax.set_ylabel('[A]')
#        ax.grid(True)
#        ax.legend(labels = ['Curvature'])
            
    return utmPoseX, utmPoseY, heading, np.asarray(Roll), np.asarray(Pitch), \
np.asarray(Yaw), np.asarray(Current), Speed, traversedDist, segmentPath, \
GPSspeed, dT, np.asarray(Time), np.asarray(Slope)


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
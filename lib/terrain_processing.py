# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 10:33:42 2019

@author: Richi
"""
import numpy as np
from numpy import array
from scipy import signal
import sys

rad2deg = 180/np.pi

# =============================================================================
#     Get Convoluted Slope Maps: Gradient Map and Aspect Maps
# =============================================================================
def getConvSlopeMaps(XX, YY, DEM, occupancyRadius, res, globalRes):
    r = int((occupancyRadius+globalRes)/res)
    r = r + 1 - r%2
    y,x = np.ogrid[-r: r+1, -r: r+1]
    convMatrix = x**2+y**2 <= r**2
    convMatrix = convMatrix.astype(float)
    convMatrix = convMatrix/convMatrix.sum()

    smoothDEM = signal.medfilt2d(DEM,r)
    smoothDEM = signal.convolve2d(DEM, convMatrix, mode='same',\
                                  boundary='symm')
    smoothDEM = signal.convolve2d(smoothDEM, convMatrix, mode='same',\
                                  boundary='symm')
    nX, nY, nZ = surface_normal(XX,YY,smoothDEM)
    gradientMap = rad2deg*np.arccos(nZ)
#    gradientMap = np.maximum(gradientMap,signal.convolve2d(gradientMap, convMatrix, mode='same',\
#                                  boundary='symm'))
    aspectMapX = nX/np.sqrt(nX**2+nY**2)
    aspectMapY = nY/np.sqrt(nX**2+nY**2)
#    dX,dY = np.gradient(smoothDEM,*[res, res])
#    gradientMap = np.zeros_like(dX)
#    gradientMap = rad2deg*np.arctan(np.sqrt(dX**2+dY**2))
#    dnX = dX/np.sqrt(dX**2+dY**2)
#    dnY = dY/np.sqrt(dX**2+dY**2)
#    aspectMap = -rad2deg*np.arctan2(dnX,-dnY)
    
    return gradientMap, aspectMapX, aspectMapY, smoothDEM

def surface_normal(x,y,z):

# =============================================================================
#    Obtain surface normals in every point of a meshgrid
# =============================================================================

    xx = x
    yy = y
    zz = z
    m,n = x.shape
    m = int(m)
    n = int(n)

    stencil1 = array([[0,0,0],[1,0,-1],[0,0,0]])/2
    stencil2 = array([[0,-1,0],[0,0,0],[0,1,0]])/2  
    
    
    
    xx = np.vstack((3*xx[0,:]-3*xx[1,:]+xx[2,:],xx,3*xx[m-1,:]-3*xx[m-2,:]+xx[m-3,:]))
    xx = np.hstack((array([3*xx[:,0]-3*xx[:,1]+xx[:,2]]).T,xx,array([3*xx[:,n-1]-3*xx[:,n-2]+xx[:,n-3]]).T))
    yy = np.vstack((3*yy[0,:]-3*yy[1,:]+yy[2,:],yy,3*yy[m-1,:]-3*yy[m-2,:]+yy[m-3,:]))
    yy = np.hstack((array([3*yy[:,0]-3*yy[:,1]+yy[:,2]]).T,yy,array([3*yy[:,n-1]-3*yy[:,n-2]+yy[:,n-3]]).T))
    zz = np.vstack((3*zz[0,:]-3*zz[1,:]+zz[2,:],zz,3*zz[m-1,:]-3*zz[m-2,:]+zz[m-3,:]))
    zz = np.hstack((array([3*zz[:,0]-3*zz[:,1]+zz[:,2]]).T,zz,array([3*zz[:,n-1]-3*zz[:,n-2]+zz[:,n-3]]).T))
    
    
    ax = -signal.convolve2d(xx, np.flip(stencil1,0), mode='valid')
    ay = -signal.convolve2d(yy, np.flip(stencil1,0), mode='valid')
    az = -signal.convolve2d(zz, np.flip(stencil1,0), mode='valid')

    bx = signal.convolve2d(xx, np.flip(stencil2,0), mode='valid')
    by = signal.convolve2d(yy, np.flip(stencil2,0), mode='valid')
    bz = signal.convolve2d(zz, np.flip(stencil2,0), mode='valid')
    
    nx = -(ay*bz - az*by)
    ny = -(az*bx - ax*bz)
    nz = -(ax*by - ay*bx)
    
    mag = np.sqrt(nx*nx+ny*ny+nz*nz)
    mag[np.where(mag == 0)] = sys.float_info.epsilon
    nxout = nx/mag
    nyout = ny/mag
    nzout = nz/mag
    
    return nxout,nyout,nzout
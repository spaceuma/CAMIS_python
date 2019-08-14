# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 12:30:57 2019

@author: Richi
"""

import math
import numpy as np


def getHexGrid(Xmap,Ymap,res):
    DX = Xmap[-1,-1] - Xmap[0,0]
    DY = Ymap[-1,-1] - Ymap[0,0]
    JMax = math.ceil(2*DY/(math.sqrt(3)*res))
    IMax = math.ceil((DX+DY/math.sqrt(3))/res)
    II,JJ = np.meshgrid(np.linspace(0,IMax,IMax+1),np.linspace(0,JMax,JMax+1))
    hX = Xmap[0,0] + res*(II + .5*JJ) - DY/math.sqrt(3)
    hY = Ymap[0,0] + res*math.sqrt(3)/2*JJ
    XX,YY = np.meshgrid(np.linspace(0,np.ceil(DX),np.ceil(DX)+1),np.linspace(0,np.ceil(DY),np.ceil(DY)+1))
    xy2J = 2*YY/(np.sqrt(3)*res)
    xy2I = (DY/np.sqrt(3)+ XX)/res-0.5*xy2J
    return hX, hY, xy2I, xy2J
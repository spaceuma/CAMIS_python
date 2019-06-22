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
    return hX, hY
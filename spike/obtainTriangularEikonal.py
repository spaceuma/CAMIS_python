# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:15:02 2020

@author: rsanchez
"""

import numpy as np

h = 1.0

P = np.array([[ .5, np.sqrt(3)/2],
              [-.5, np.sqrt(3)/2]])

a = np.array([[1/h],[1/h]])

Q = np.linalg.inv(P.dot(P.transpose()))

M = ((a.transpose()).dot(Q)).dot(a)

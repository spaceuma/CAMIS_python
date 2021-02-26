# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 10:16:31 2020

@author: rsanchez
"""

import numpy as np
import matplotlib.pyplot as plt
from context import camis
import copy
import yaml
import matplotlib.gridspec as gridspec
try:
    from scipy import signal
except:
    raise ImportError('ERROR: scipy module could not be imported')

with open("spikeCAMIS.yml", 'r') as file:
    robot1 = yaml.full_load(file) 
r1 = camis.CamisDrivingModel(robot1) 
r1.showBraking()
r1.showDirCosts()
r1.showAnisotropy()
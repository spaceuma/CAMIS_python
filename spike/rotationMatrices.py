# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 18:09:28 2020

@author: Richi
"""

from scipy.spatial.transform import Rotation as R

r = R.from_euler('zyx', [-55.1, -15.2, -5.6], degrees=True)
r3 = R.from_euler('zyx', [00.0, 0.0, 90.0], degrees=True)
q = r.inv()

r2 = R.from_euler('zyx', [38.8, -9.8, -10.4 ], degrees=True)

p = q*r2
yaw,pitch,roll = p.as_euler('zyx', degrees=True)
print(roll)
print(pitch)
print(yaw)

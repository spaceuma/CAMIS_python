# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import optimize
from scipy.optimize import least_squares

deg2rad = np.pi/180
rad2deg = 180/np.pi

steepnessArray = np.linspace(0,90,200)

slip_ratioA = 0.07 * np.exp(0.10 * steepnessArray)
slip_angleA = 1.32 * np.exp(0.16 * steepnessArray)
slip_ratioB = 0.04 * np.exp(0.07 * steepnessArray)
slip_angleB = 1.20 * np.exp(0.08 * steepnessArray)


plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, ax = plt.subplots(figsize=(3, 3),constrained_layout=True)
p1, = ax.plot(steepnessArray, slip_ratioA)
p3, = ax.plot(steepnessArray, slip_ratioB)
p5, = ax.plot([0], marker='None',
           linestyle='None', label='Model A')
model_labels = ['Wheel Model', 'Track Model']
ax.set_ylim([0.0, 1.0])
ax.set_xlim([0.0, 90.0])
ax.set_xlabel('Steepness α [degrees]')
ax.set_ylabel('Slip Ratio $s_r$')
#l1 = ax.legend((p1[0],p2[0]),('$s_r = 0.07 e^{0.10 α}$', \
#          '$\cos_{s_a} = \cos (1.32 e^{0.16 α} π / 180.0)$'), loc = 2)
#l2 = ax.legend((p3[0],p4[0]),('$s_r = 1.32 e^{0.16 α}$', \
#          '$\cos_{s_a} = \cos (0.20 e^{0.08 α} π / 180.0)$'), loc = 7)
#plt.gca().add_artist(l1)
l1 = ax.legend([p5,p1,p5,p5,p3], [r'$Wheel \ Model$'] + \
               ['$s_r = 0.07 e^{0.10 α}$'] + [''] + \
                [r'$Track \ Model$'] + \
                ['$s_r = 0.04 e^{0.07 α}$'])
plt.minorticks_on()
plt.grid(b=True,which='minor', linestyle = '--')
plt.grid(b=True,which='major', linewidth = 1)
plt.show()

fig, ax = plt.subplots(figsize=(3, 3),constrained_layout=True)
p2, = ax.plot(steepnessArray, slip_angleA)
p4, = ax.plot(steepnessArray, slip_angleB)
p5, = ax.plot([0], marker='None',
           linestyle='None', label='Wheel Model')
l1 = ax.legend([p5,p2,p5,p5,p4], [r'$Wheel \ Model$'] + \
               ['$s_a = 1.32 e^{0.16 α}$'] + [''] + \
                [r'$Track \ Model$'] + \
                ['$s_a = 1.20 e^{0.08 α}$'])
ax.set_ylim([0.0, 90.0])
ax.set_xlim([0.0, 90.0])
ax.set_xlabel('Steepness α [degrees]')
ax.set_ylabel('Slip Angle $s_a$ [degrees]')
plt.minorticks_on()
plt.grid(b=True,which='minor', linestyle = '--')
plt.grid(b=True,which='major', linewidth = 1)
plt.show()

fig, ax = plt.subplots(figsize=(5, 3),constrained_layout=True)
ax.plot(steepnessArray, 1/(1-slip_ratioA))
ax.plot(steepnessArray, 1/(1-slip_ratioB))
ax.plot(steepnessArray, 1/np.cos(slip_angleA*deg2rad))
ax.plot(steepnessArray, 1/np.cos(slip_angleB*deg2rad))
ax.set_ylim([0.0, 5.0])
ax.set_xlabel('Steepness [degrees]')
ax.set_ylabel('Value')
ax.legend(('Slip Ratio Model A', 'Slip Ratio Model B'))
plt.show()

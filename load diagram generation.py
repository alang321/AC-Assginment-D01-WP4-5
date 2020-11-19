#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 15:51:54 2020

@author: olivierheu
"""


import numpy as np 
import matplotlib.pyplot as plt
from scipy import interpolate


#test speeds (inputs)

V_A = 230
V_D = 300
V_F = 270
V_S1 = 200

#test loads (inputs)
n_max = 3
n_min = -1.2


speeds = [0,V_A, V_D, V_D, V_F, V_S1]  
n_values = [0,n_max, n_max, 0, n_min, n_min]

diagram_values = np.column_stack((speeds, n_values))

maneuvre_envelope = plt.figure()

ax = maneuvre_envelope.add_subplot(1, 1, 1)
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')


for i in range(1, len(speeds)):
    plt.plot(speeds[i:i+2], n_values[i:i+2], "black")
    
    
    
plt.show()






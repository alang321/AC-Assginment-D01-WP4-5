#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 15:51:54 2020

@author: olivierheu
"""


import numpy as np 
import matplotlib.pyplot as plt
from scipy import interpolate



speeds = [150, 230, 300, 300, 270, 210, 150]   # origin, V_A, V_D, V_D, V_F, V_S1, 
n_values = [0, 3, 3, 0, -1.2, -1.2, 0]


diagram_values = np.column_stack((speeds, n_values))



for i in range(1, len(speeds)):
    plt.plot(speeds[i:i+2], n_values[i:i+2], "ro-")
    
plt.show()



x = [speeds[0],speeds[5],speeds[1]]             
y=[0,2,n_values[1]]

x2 = np.linspace(x[0], x[-1], 100)
y2 = interpolate.pchip_interpolate(x, y, x2)
plt.plot(x2, y2, "ro-")
plt.plot(x, y, "o")


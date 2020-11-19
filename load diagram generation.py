#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 15:51:54 2020

@author: olivierheu
"""


import numpy as np 
import matplotlib.pyplot as plt
import math

###---MANEUVRE LOAD DIAGRAM---###

#test speeds (inputs)

V_A = 230
V_D = 300
V_F = 270
V_S1 = 200
V_S0 = 150

#test loads (inputs)
n_max = 3
n_min = -1.2

def get_V_n_2(Vs1):

    V_n_2 = Vs1 * (math.sqrt(2))

    return V_n_2


speeds = [V_A, V_D, V_D, V_F, V_S1]  
n_values = [n_max, n_max, 0, n_min, n_min]

plt.plot(speeds, n_values, 'black')
plt.title('Maneuvre Envelope')
plt.xlabel('Velocity')
plt.ylabel('Load')
plt.text(V_A, n_max, 'V_A') #add text to diagram


def f(x):
    return (x / V_S1)**2


x1 = np.linspace(0,V_A, 1000)



plt.plot(x1, (x1 / V_S0)**2, 'black')  # (0,0) to V_A curve
#plt.plot(x1, min((x1 / V_S0)**2, 2), 'black')  #flaps down curve
 
plt.show()






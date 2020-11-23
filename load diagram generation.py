#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 15:51:54 2020

@author: olivierheu
"""


import numpy as np 
import matplotlib.pyplot as plt
import math
import loadingdiagrams

###---MANEUVRE LOAD DIAGRAM---###
def plot_maneuver(V_A, V_D, V_F, V_S0, V_S1):

    def f(x):
        return (x / V_S1)**2
                  
    speeds = [V_A, V_D, V_D, V_F, V_S1]  
    n_values = [n_max, n_max, 0, n_min, n_min]

    plt.plot(speeds, n_values, 'black')
    plt.title('Maneuvre Envelope')
    plt.xlabel('Velocity')
    plt.ylabel('Load')
    plt.text(V_A, n_max, 'V_A') #add text to diagram


    x1 = np.linspace(0,V_A, 1000)
    x2 = np.linspace(0, V_S1 * math.sqrt(2), 1000)
    x3 = np.linspace(0, V_S1, 1000)
    y2 = []      #flaps down curve n values

    for i in x2:
        a = (i / V_S0)**2
        a = min(a, 2)
        y2.append(a)

    plt.plot(x1, f(x1), 'black')  # (0,0) to V_A curve
    plt.plot(x2, y2, 'black')  #flaps down curve
    plt.plot(x3, -f(x3), 'black')
 
    return plt.show()

for sets in V_all_list_sorted:

    V_A = sets[0]
    V_D = sets[2]
    V_F = sets[3]
    V_S0 = sets[4]
    V_S1 = sets[5]

    plot_maneuver(V_A, V_D, V_F, V_S0, V_S1)





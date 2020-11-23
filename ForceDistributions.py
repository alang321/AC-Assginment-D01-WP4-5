import AerodynamicLoading
from aircraftProperties import AircraftProperties
import numpy as np
import scipy as sp
from scipy import integrate
from math import sin, cos, radians, degrees
import matplotlib.pyplot as plt
from scipy import interpolate

def a(x):
    return AerodynamicLoading.Liftacc(x)
def b (x):
    return AerodynamicLoading.Dragacc(x)
def c (x):
    return AerodynamicLoading.Momentacc(x)
def d (x):
    return AerodynamicLoading.Liftacc10(x)
def e (x):
    return AerodynamicLoading.Dragacc10(x)
def f (x):
    return AerodynamicLoading.Momentacc10(x)


Lift = integrate.quad(a, 0, AircraftProperties.Planform["span"]/2)
Drag = integrate.quad(b, 0, AircraftProperties.Planform["span"]/2)
Moment = integrate.quad(c, 0, AircraftProperties.Planform["span"]/2)
Lift10 = integrate.quad(d, 0, AircraftProperties.Planform["span"]/2)
Drag10 = integrate.quad(e, 0, AircraftProperties.Planform["span"]/2)
Moment10 = integrate.quad(f, 0, AircraftProperties.Planform["span"]/2)

AoA0 = 0      # deg
AoA10 = 10
NormalF0lst = []
NormalF10lst = []
TangentialF0lst = []
TangentialF10lst = []
for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
    NormalF0lst.append(cos(radians(AoA0))*AerodynamicLoading.Liftacc(i) + sin(radians(AoA0))*AerodynamicLoading.Dragacc(i))
    NormalF10lst.append(cos(radians(AoA10))*AerodynamicLoading.Liftacc10(i) + sin(radians(AoA10))*AerodynamicLoading.Dragacc10(i))
    TangentialF0lst.append(-sin(radians(AoA0))*AerodynamicLoading.Liftacc(i) + cos(radians(AoA0))*AerodynamicLoading.Dragacc(i))
    TangentialF10lst.append(-sin(radians(AoA10))*AerodynamicLoading.Liftacc10(i) + cos(radians(AoA10))*AerodynamicLoading.Dragacc10(i))

NormalF0 = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), NormalF0lst, kind='cubic', fill_value="extrapolate")
NormalF10 = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), NormalF10lst, kind='cubic', fill_value="extrapolate")
TangentialF0 = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), TangentialF0lst, kind='cubic', fill_value="extrapolate")
TangentialF10 = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), TangentialF10lst, kind='cubic', fill_value="extrapolate")

print(AerodynamicLoading.Liftacc(3), AerodynamicLoading.Liftacc10(3), AerodynamicLoading.Dragacc(3), AerodynamicLoading.Dragacc10(3))
print(NormalF0(3), NormalF10(3), TangentialF0(3), TangentialF10(3))

def g (x):
    return NormalF0(x)
def h (x):
    return NormalF10(x)

ShearF0lst = []
for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
    intShear = integrate.quad(g, i, AircraftProperties.Planform["span"]/2)
    if i >= 10.02:    
        ShearF0lst.append(intShear[0])
    else: 
        ShearF0lst.append(intShear[0] - AircraftProperties.Engine["weight"])
ShearF0 = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), ShearF0lst, kind='cubic', fill_value="extrapolate")

ShearF10lst = []
for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
    intShear10 = integrate.quad(h, i, AircraftProperties.Planform["span"]/2)
    if i >= 10.02:    
        ShearF10lst.append(intShear10[0])
    else: 
        ShearF10lst.append(intShear10[0] - AircraftProperties.Engine["weight"])
ShearF10 = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), ShearF10lst, kind='cubic', fill_value="extrapolate")


def drawgraphs():
    plt.subplot(221)
    X = []
    Y0 = []
    Y10 = []
    for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
        X.append(i)
        Y0.append(ShearF0(i))
        Y10.append(ShearF10(i))
    plt.plot(X, Y0, color="red")
    plt.plot(X, Y10, color="blue")
    plt.ylabel("shear")
    plt.show()


drawgraphs()

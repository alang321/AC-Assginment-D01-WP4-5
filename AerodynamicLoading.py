import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
from aircraftProperties import AircraftProperties
import math

rawdata0 = np.genfromtxt("data/MainWing_a=0.00_v=10.00ms.txt", skip_header=40, max_rows=19)
rawdata0 = rawdata0.T
rawdata10 = np.genfromtxt("data/MainWing_a=10.00_v=10.00ms.txt", skip_header=40, max_rows=19)
rawdata10 = rawdata10.T
# print(rawdata0)
# print(rawdata10)

# Makes all data into lists
Ylst = rawdata0[0]
Chordlst = rawdata0[1]
Cl0lst = rawdata0[3]
Cl10lst = rawdata10[3]
ICd0lst = rawdata0[5]
ICd10lst = rawdata10[5]
Cm0lst = rawdata0[7]
Cm10lst = rawdata10[7]

# Makes lists into functions
Chord = sp.interpolate.interp1d(Ylst, Chordlst, kind='cubic', fill_value="extrapolate")
Cl0 = sp.interpolate.interp1d(Ylst, Cl0lst, kind='cubic', fill_value="extrapolate")
Cl10 = sp.interpolate.interp1d(Ylst, Cl10lst, kind='cubic', fill_value="extrapolate")
ICd0 = sp.interpolate.interp1d(Ylst, ICd0lst, kind='cubic', fill_value="extrapolate")
ICd10 = sp.interpolate.interp1d(Ylst, ICd10lst, kind='cubic', fill_value="extrapolate")
Cm0 = sp.interpolate.interp1d(Ylst, Cm0lst, kind='cubic', fill_value="extrapolate")
Cm10 = sp.interpolate.interp1d(Ylst, Cm10lst, kind='cubic', fill_value="extrapolate")

rho = 1.225
V = 10
q = 1 / 2 * rho * (V ** 2)
Liftacclst = []
for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
    Liftacclst.append(Cl0(i) * Chord(i) * q)
Liftacc = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), Liftacclst, kind='cubic', fill_value="extrapolate")
Dragacclst = []
for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
    Dragacclst.append(ICd0(i) * Chord(i) * q)
Dragacc = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), Dragacclst, kind='cubic', fill_value="extrapolate")
Momentacclst = []
for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
    Momentacclst.append(Cm0(i) * Chord(i) * q)
Momentacc = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), Momentacclst, kind='cubic', fill_value="extrapolate")

# 10 degrees AoA
Liftacc10lst = []
for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
    Liftacc10lst.append(Cl10(i) * Chord(i) * q)
Liftacc10 = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), Liftacc10lst, kind='cubic', fill_value="extrapolate")
Dragacc10lst = []
for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
    Dragacc10lst.append(ICd10(i) * Chord(i) * q)
Dragacc10 = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), Dragacc10lst, kind='cubic', fill_value="extrapolate")
Momentacc10lst = []
for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
    Momentacc10lst.append(Cm10(i) * Chord(i) * q)
Momentacc10 = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), Momentacc10lst, kind='cubic',
                                      fill_value="extrapolate")

# Distribution at random AoA / Coeff
CL_0 = 0.430299
CL_10 = 1.252323
CL_d = 2  # Desired distribution for this coefficcient
CL_dList = []
for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
    CL_dList.append(Liftacc(i) + ((CL_d - CL_0) / (CL_10 - CL_0)) * (Liftacc10(i) - Liftacc(i)))
CL_dacc = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), CL_dList, kind='cubic', fill_value="extrapolate")

Cd_0 = 0.006490
Cd_10 = 0.054091
Cd_d = 0.02  # Desired distribution for this coefficcient
Cd_dList = []
for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
    Cd_dList.append(Dragacc(i) + ((Cd_d - Cd_0) / (Cd_10 - Cd_0)) * (Dragacc10(i) - Dragacc(i)))
Cd_dacc = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), Cd_dList, kind='cubic', fill_value="extrapolate")

Cm_0 = -0.621253
Cm_10 = -1.60557
Cm_d = -1  # Desired distribution for this coefficcient
Cm_dList = []
for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
    Cm_dList.append(Momentacc(i) + ((Cm_d - Cm_0) / (Cm_10 - Cm_0)) * (Momentacc10(i) - Momentacc(i)))
Cm_dacc = sp.interpolate.interp1d(np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1), Cm_dList, kind='cubic', fill_value="extrapolate")

# calculation AoA
AoAL = math.degrees(math.asin(((CL_d - CL_0) / (CL_10 - CL_0)) * math.sin(math.radians(10))))
AoAD = math.degrees(math.asin(((Cd_d - Cd_0) / (Cd_10 - Cd_0)) * math.sin(math.radians(10))))
AoAM = math.degrees(math.asin(((Cm_d - Cm_0) / (Cm_10 - Cm_0)) * math.sin(math.radians(10))))

# print(AoAL)
# print(AoAD)
# print(AoAM)


def drawgraphs():
    plt.subplot(221)
    X = []
    Y0 = []
    Y10 = []
    Yd = []
    for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
        X.append(i)
        Y0.append(Liftacc(i))
        Y10.append(Liftacc10(i))
        Yd.append(CL_dacc(i))
    plt.plot(X, Y0, color="red")
    plt.plot(X, Y10, color="blue")
    plt.plot(X, Yd, color="green")
    plt.ylabel("Lift")

    plt.subplot(222)
    X = []
    Y0 = []
    Y10 = []
    Yd = []
    for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
        X.append(i)
        Y0.append(Dragacc(i))
        Y10.append(Dragacc10(i))
        Yd.append(Cd_dacc(i))
    plt.plot(X, Y0, color="red")
    plt.plot(X, Y10, color="blue")
    plt.plot(X, Yd, color="green")
    plt.xlabel("Span Position")
    plt.ylabel("Drag")

    plt.subplot(223)
    X = []
    Y0 = []
    Y10 = []
    Yd = []
    for i in np.arange(0, AircraftProperties.Planform["span"] / 2, 0.1):
        X.append(i)
        Y0.append(Momentacc(i))
        Y10.append(Momentacc10(i))
        Yd.append(Cm_dacc(i))
    plt.plot(X, Y0, color="red")
    plt.plot(X, Y10, color="blue")
    plt.plot(X, Yd, color="green")
    plt.xlabel("Span Position")
    plt.ylabel("Moment")

    plt.show()
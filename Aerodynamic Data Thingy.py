import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
from aircraftProperties import AircraftProperties

rawdata0 = np.genfromtxt("data/MainWing_a=0.00_v=10.00ms.txt", skip_header=40, max_rows=19)
rawdata0 = rawdata0.T
rawdata10 = np.genfromtxt("data/MainWing_a=10.00_v=10.00ms.txt", skip_header=40, max_rows=19)
rawdata10 = rawdata10.T
# print(rawdata0)
# print(rawdata10)

# Makes all data into lists
Ylst = rawdata0[0]
Cl0lst = rawdata0[3]
Cl10lst = rawdata10[3]
ICd0lst = rawdata0[5]
ICd10lst = rawdata10[5]
Cm0lst = rawdata0[7]
Cm10lst = rawdata10[7]

# Makes lists into functions
Cl0 = sp.interpolate.interp1d(Ylst, Cl0lst, kind='cubic', fill_value="extrapolate")
Cl10 = sp.interpolate.interp1d(Ylst, Cl10lst, kind='cubic', fill_value="extrapolate")
ICd0 = sp.interpolate.interp1d(Ylst, ICd0lst, kind='cubic', fill_value="extrapolate")
ICd10 = sp.interpolate.interp1d(Ylst, ICd10lst, kind='cubic', fill_value="extrapolate")
Cm0 = sp.interpolate.interp1d(Ylst, Cm0lst, kind='cubic', fill_value="extrapolate")
Cm10 = sp.interpolate.interp1d(Ylst, Cm10lst, kind='cubic', fill_value="extrapolate")

X = []
Y = []
for i in np.arange(0, AircraftProperties.Planform["span"]/2, 0.1):
    X.append(i)
    Y.append(Cl0(i))
plt.plot(X, Y, color="red")
plt.ylim(0)
plt.show()
import numpy as np
import scipy as sp
from scipy import interpolate


rawdatashear = np.genfromtxt("data/bucklingcoeff.dat", max_rows=7)
rawdatashear = rawdatashear.T


rawdatacompression = np.genfromtxt("data/bucklingcoeff.dat",skip_header=9, max_rows=20)
rawdatacompression = rawdatacompression.T

abshear = rawdatashear[0]
abcompression = rawdatacompression[0]
coeffshear = rawdatashear[1]
coeffcompression = rawdatacompression[1]

shearfunc = sp.interpolate.interp1d(abshear, coeffshear, kind='linear', fill_value="extrapolate")
compressionfunc = sp.interpolate.interp1d(abcompression, coeffcompression, kind='linear', fill_value="extrapolate")


#print(shearfunc(1.7))
#print(compressionfunc(1.3))


def getbucklingcoeffshear(ab):
    return shearfunc(ab)

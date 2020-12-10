import numpy as np
import scipy as sp
from scipy import interpolate


rawdatashearhinge = np.genfromtxt("data/bucklingcoeff.dat", max_rows=7)
rawdatashearhinge = rawdatashearhinge.T


rawdatashearclamp = np.genfromtxt("data/bucklingcoeff.dat",skip_header=29,  max_rows=9)
rawdatashearclamp = rawdatashearclamp.T

rawdatacompressionSS = np.genfromtxt("data/bucklingcoeff.dat",skip_header=9, max_rows=20)
rawdatacompressionSS = rawdatacompressionSS.T

rawdatacompressionC = np.genfromtxt("data/bucklingcoeff.dat",skip_header=39, max_rows=14)
rawdatacompressionC = rawdatacompressionC.T

rawdatacompressionCC = np.genfromtxt("data/bucklingcoeff.dat",skip_header=54, max_rows=19)
rawdatacompressionCC = rawdatacompressionCC.T

abshearH = rawdatashearhinge[0]
abcompressionSS = rawdatacompressionSS[0]
coeffshearH = rawdatashearhinge[1]
coeffcompressionSS = rawdatacompressionSS[1]
abshearC = rawdatashearclamp[0]
coeffshearC = rawdatashearclamp[1]
abcompressionC = rawdatacompressionC[0]
coeffcompressionC = rawdatacompressionC[1]
coeffcompressionCC = rawdatacompressionCC[1]
abcompressionCC = rawdatacompressionCC[0]

shearfuncH = sp.interpolate.interp1d(abshearH, coeffshearH, kind='linear', fill_value="extrapolate")
compressionfuncH = sp.interpolate.interp1d(abcompressionSS, coeffcompressionSS, kind='linear', fill_value="extrapolate")
shearfuncC = sp.interpolate.interp1d(abshearC, coeffshearC, kind='linear', fill_value="extrapolate")
compressionfuncC = sp.interpolate.interp1d(abcompressionC, coeffcompressionC, kind='linear', fill_value="extrapolate")
compressionfuncCC = sp.interpolate.interp1d(abcompressionCC, coeffcompressionCC, kind='linear', fill_value="extrapolate")


print(compressionfuncCC(4.9))



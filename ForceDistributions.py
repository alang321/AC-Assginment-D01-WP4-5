import AerodynamicLoading
import numpy as np
import scipy as sp
from scipy import integrate

def a (x):
    return AerodynamicLoading.Liftacc(x)


y = integrate.quad(a, 0, 8)

print(y)

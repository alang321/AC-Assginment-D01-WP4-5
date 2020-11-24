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

# wing coefficients for 0 at 10 deg aoa
CL_0 = 0.430299
CL_10 = 1.252323
Cd_0 = 0.006490
Cd_10 = 0.054091
Cm_0 = -0.621253
Cm_10 = -1.60557

# Makes lists into functions
# coefficeint functions at span location mach <0.2
Chord = sp.interpolate.interp1d(Ylst, Chordlst, kind='cubic', fill_value="extrapolate")
Cl0Func = sp.interpolate.interp1d(Ylst, Cl0lst, kind='cubic', fill_value="extrapolate")
Cl10Func = sp.interpolate.interp1d(Ylst, Cl10lst, kind='cubic', fill_value="extrapolate")
Cd0Func = sp.interpolate.interp1d(Ylst, ICd0lst, kind='cubic', fill_value="extrapolate")
Cd10Func = sp.interpolate.interp1d(Ylst, ICd10lst, kind='cubic', fill_value="extrapolate")
Cm0Func = sp.interpolate.interp1d(Ylst, Cm0lst, kind='cubic', fill_value="extrapolate")
Cm10Func = sp.interpolate.interp1d(Ylst, Cm10lst, kind='cubic', fill_value="extrapolate")

coefficientFunctions = [[Cl0Func, Cl10Func],
                        [Cd0Func, Cd10Func],
                        [Cm0Func, Cm10Func]]

#for drawing
aerodynamicCoefficientsDrawing = [[Cl0Func, "Lift coefficient at 0 deg"], [Cl10Func, "Lift coefficient at 10 deg"], [Cd0Func, "Drag coefficient at 0 deg"], [Cd10Func, "Drag coefficient at 10 deg"], [Cm0Func, "Moment around c/4 0 deg"], [Cm10Func, "Moment around c/4 10 deg"]]


def getNormalTangentialMomentAOA(cLd, xovercCentroid, v, altitude):
    temp, p, rho = getISAParameters(altitude)
    machNumber = v/(1.4*287.05*temp)**0.5
    q = 1 / 2 * rho * (v**2)

    if machNumber > 0.2:
        adjustedCoefficients = getCompressionAdjustedCoefficients(machNumber)
    else:
        adjustedCoefficients = [[CL_0, CL_10], [Cd_0, Cd_10], [Cm_0, Cm_10]]

    aoa = getAOAatCL(cLd, adjustedCoefficients)

    CDd = getCoefficientatAOA(aoa, 1, adjustedCoefficients)
    CMd = getCoefficientatAOA(aoa, 2, adjustedCoefficients)

    clFreestream = getCoefficientDistribution(cLd, 0, machNumber, adjustedCoefficients)
    cdFreestream = getCoefficientDistribution(CDd, 1, machNumber, adjustedCoefficients)
    cm = getCoefficientDistribution(CMd, 2, machNumber, adjustedCoefficients)
    # for shit of moment equation see first line of https://brightspace.tudelft.nl/d2l/le/content/292972/viewContent/1906402/View
    cmx = lambda y: cm(y) + (xovercCentroid - 0.25) * clFreestream(y)

    liftFreestream = lambda y: clFreestream(y) * q * Chord(y)
    dragFreestream = lambda y: cdFreestream(y) * q * Chord(y)

    Moment = lambda y: cmx(y) * q * Chord(y)

    Lift = lambda y: np.cos(aoa) * liftFreestream(y) + np.sin(aoa) * dragFreestream(y)
    Drag = lambda y: -liftFreestream(y) * np.sin(aoa) + np.cos(aoa) * dragFreestream(y)

    return Lift, Drag, Moment, aoa

def getCoefficientDistribution(desiredCoefficient, index, machNumber, adjustedCoefficients):
    coeff0 = adjustedCoefficients[index][0]
    coeff10 = adjustedCoefficients[index][1]
    func0 = coefficientFunctions[index][0]
    func10 = coefficientFunctions[index][1]

    if machNumber > 0.2:
        return lambda y: (func0(y) + ((desiredCoefficient - coeff0) / (coeff10 - coeff0)) * (func10(y) - func0(y))) * (1 / getBeta(machNumber))
    else:
        return lambda y: (func0(y) + ((desiredCoefficient - coeff0) / (coeff10 - coeff0)) * (func10(y) - func0(y)))


def getCoefficientatAOA(aoa, index, adjustedCoefficients):
    return ((np.sin(aoa)*(adjustedCoefficients[index][1]-adjustedCoefficients[index][0]))/np.sin(np.deg2rad(10)))+adjustedCoefficients[index][0]


# calculation AoA
def getAOAatCL(CL_d, adjustedCoefficients):
    return np.arcsin(((CL_d - adjustedCoefficients[0][0]) / (adjustedCoefficients[0][1] - adjustedCoefficients[0][0])) * np.sin(np.deg2rad(10)))


def getBeta(machNumber):
    return (abs(1-machNumber**2))**0.5


def getCompressionAdjustedCoefficients(machNumber):
    cd0 = AircraftProperties.Planform["cd0 wing"]
    zeroliftAOA = np.deg2rad(AircraftProperties.Airfoil["zero lift aoa"])
    e = AircraftProperties.Planform["e"]
    A = AircraftProperties.Planform["aspect ratio"]

    beta = getBeta(machNumber)
    eff = 0.95
    sweep = np.deg2rad(AircraftProperties.Planform["half-chord sweep"])

    clalpha = (2 * np.pi * A) / (2 + np.sqrt(4 + (((A * beta) ** 2) / eff) * (1 + ((np.tan(sweep) ** 2) / (beta ** 2)))))

    d = clalpha * -zeroliftAOA

    CL0adjusted = d
    CL10adjusted = clalpha * np.deg2rad(10) + d

    Cm0adjusted = Cm_0 * (CL0adjusted/CL_0)
    Cm10adjusted = Cm_10 * (CL10adjusted/CL_10)

    Cd0adjusted = cd0 + (CL0adjusted**2)/(np.pi * A * e)
    Cd10adjusted = cd0 + (CL10adjusted**2)/(np.pi * A * e)

    return [[CL0adjusted, CL10adjusted], [Cd0adjusted, Cd10adjusted], [Cm0adjusted, Cm10adjusted]]


def getISAParameters(h):
    g = 9.80665
    R = 287.05
    T0 = 288.15
    p0 = 101325.0
    a = -0.0065

    T1 = T0 + (a * h)
    p1 = p0 * ((T1 / T0) ** (-(g / (a * R))))
    rho = p1 / (R * T1)

    return T1, p1, rho

def drawAerodynamicCoefficients(fidelity=50):
    a = np.linspace(0, AircraftProperties.Planform["span"]/2, fidelity)
    fig, plots = plt.subplots(3, 2)
    fig.suptitle('Aerodynamic Coefficients\n', fontsize=16)
    for row in range(3):
        for col in range(2):
            plots[row, col].plot(a, [aerodynamicCoefficientsDrawing[row * 2 + col][0](i) for i in a])
            plots[row, col].set_title(aerodynamicCoefficientsDrawing[row * 2 + col][1])
            plots[row, col].set(xlabel='semi-span [m]')
    fig.tight_layout(pad=0.2)

    plt.show()

def drawAerodynamicForces(aerodynamicForces, fidelity=50):
    a = np.linspace(0, AircraftProperties.Planform["span"]/2, fidelity)
    plotNames = ["Normal", "Tangential", "Distributed Moment at x/c"]
    plotUnits = ["N", "N", "Nm"]
    fig, plots = plt.subplots(3)
    fig.suptitle('Aerodynamic Forces\n', fontsize=16)
    for row in range(3):
        plots[row].plot(a, [aerodynamicForces[row](i) for i in a])
        plots[row].set_title(plotNames[row])
        plots[row].set(xlabel='semi-span [m]', ylabel=plotUnits[row])

    fig.tight_layout(pad=0.2)

    plt.show()

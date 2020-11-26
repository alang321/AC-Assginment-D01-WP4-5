from aircraftProperties import AircraftProperties
from Polygon import StringerType
from wingbox import Wingbox
import AerodynamicLoading
import numpy as np
import scipy as sp
from scipy import interpolate
from WingLoads import WingLoads
from WingLoads import PointForce


def getWingLoading(velocity, altitude, weight, loadFactor, engineThrustFactor, fuelMassDistribution=lambda y: 0, structuralMassDistribution=lambda y: 0, fuelMassFactor=1, neglectTangential=True):
    temp, p, rho = AerodynamicLoading.getISAParameters(altitude)

    S = AircraftProperties.Planform["surface area"]
    cL = (weight * loadFactor) / (1 / 2 * rho * (velocity**2) * S)
    forces = AerodynamicLoading.getNormalTangentialMomentAOA(cL, 0.375, velocity, altitude)

    aoa = forces[3]

    thrust = AircraftProperties.Engine["thrust"] * 1000 * engineThrustFactor

    # 7.369915250148878 is the chord at 10.02
    engine = PointForce([-7.369915250148878 * (0.375 + 0.2), 10.02, -1.5],
                        [-thrust, 0, -AircraftProperties.Engine["weight"]*loadFactor], "Engine")
    engineRotated = engine.getRotatedVectorAroundAxis(-aoa, [0, 1, 0])

    #add fuel and structural weight
    normal = lambda y: forces[0](y) - (fuelMassDistribution(y) + fuelMassFactor + structuralMassDistribution(y)) * np.cos(aoa) * 9.80665 * loadFactor

    if neglectTangential:
        tangential = lambda y: 0
    else:
        tangential = lambda y: forces[1](y) - (fuelMassDistribution(y) + fuelMassFactor + structuralMassDistribution(y)) * np.sin(aoa) * 9.80665 * loadFactor

    wingloading = WingLoads([engineRotated], [tangential, None, normal], [None, forces[2], None])

    return wingloading, aoa

def plotWingLoading(loading):
    print("External Forces")
    for axis in range(3):
        loading.drawExternalForces(axis)

    print("External Moments")
    for axis in range(3):
        loading.drawExternalMoments(axis)

    print("Shear Forces")
    loading.drawShearForce(0)
    loading.drawShearForce(2)

    print("Normal Force")
    loading.drawNormalForce()

    print("Internal Moments")
    for axis in range(3):
        loading.drawInternalMoment(axis)


#v = 232.412
#h = AircraftProperties.Cruise_constants["cruise altitude"]
#rho = AircraftProperties.Cruise_constants["density at cruise"]
#S = AircraftProperties.Planform["surface area"]
#cL = 1866756 / ((1 / 2) * rho * v ** 2 * S)
#forces = AerodynamicLoading.getNormalTangentialMomentAOA(cL, 0.375, v, h)
#print("Aoa:", np.rad2deg(forces[3]))


#AerodynamicLoading.drawAerodynamicCoefficients()
#AerodynamicLoading.drawAerodynamicForces(forces)

#thrustFactor = 1
#aoa = forces[3]
#thrust = AircraftProperties.Engine["thrust"] * 1000 * thrustFactor

#engine = PointForce([-7.369915250148878 * (0.375 + 0.2), 10.02, -1.5], [-thrust, 0, -AircraftProperties.Engine["weight"]], "Engine")
#engineRotated = engine.getRotatedVectorAroundAxis(-aoa, [0, 1, 0])

#engineRotated.draw()

#wingloading = WingLoads([engineRotated], [None, None, forces[0]], [None, forces[2], None])






scale = 200
stringer = StringerType([[0, 0], [0, -3/scale], [17/scale, -3/scale], [17/scale, -20/scale], [20/scale, -20/scale], [20/scale, 0]], [[0, 0], [20/scale, 0]])
stringerTop = stringer
stringerBottom = stringer.getMirrorStringerX()





yLocations = [-1, 15, 30]
thicknesses = [-1, 0.05, 0.01]
sparThicknessFunc = sp.interpolate.interp1d(yLocations, thicknesses, kind='next', fill_value="extrapolate")

yLocations = [-1, 30]
thicknesses = [-1, 0.01]
flangeThicknessTopFunc = sp.interpolate.interp1d(yLocations, thicknesses, kind='next', fill_value="extrapolate")

yLocations = [-1, 30]
thicknesses = [-1, 0.05]
flangeThicknessBottomFunc = sp.interpolate.interp1d(yLocations, thicknesses, kind='next', fill_value="extrapolate")



# rib location : fraction location of semispan, spars: [[absstart, absend, fractchordloc, thickness], []], sparThickness: [[[Thicknesses of spar], endpoint], [[Thicknesses of spar2], endpoint2]], stringersTop = [[starty, endy, fraction of chord, stringer]] flangeThickness = same as spar thickness
wingbox = Wingbox(ribLocations=[], spars=[[0, 30, 0.15, sparThicknessFunc], [0, 30, 0.6, sparThicknessFunc], [0, 15, 0.375, sparThicknessFunc]], stringersTop=[[0, 35, 0.3, stringer]], stringersBottom=[[0, 10, 0.4, stringerBottom]], sparFlangeConnectionStringerShape=stringer, flangeThicknesses=[flangeThicknessTopFunc, flangeThicknessBottomFunc], crosssectionAmount=200)
wingbox.draw(top=False)


wingbox.drawStructuralMassDistribution()
wingbox.drawFuelMassDistribution()

wingloading = getWingLoading(velocity=232.412, altitude=AircraftProperties.Cruise_constants["cruise altitude"], weight=1866756, engineThrustFactor=1, loadFactor=1, fuelMassDistribution=wingbox.getFuelMassDistribution(), fuelMassFactor=0, structuralMassDistribution=wingbox.getStructuralMassDistribution(), neglectTangential=False)
#wingloading = getWingLoading(velocity=232.412, altitude=AircraftProperties.Cruise_constants["cruise altitude"], weight=1866756, engineThrustFactor=1, loadFactor=1, neglectTangential=False)
plotWingLoading(wingloading[0])


wingbox.drawInertias()


posY = 5
wingbox.drawCrosssection(posY, drawCentroid=True, drawSidewallCenterlines=True, drawNormalStress=False)


max, min = wingbox.getMaximumMinimumNormalStress()
print("Max stress", max[0], "at y:", max[1].yLocation, "at point:", max[2])
wingbox.drawCrosssection(max[1].yLocation, drawNormalStress=True)
print("Min stress", min[0], "at y:", min[1].yLocation, "at point:", min[2])
wingbox.drawCrosssection(min[1].yLocation, drawNormalStress=True)


#deflection
wingbox.drawDeflection(wingbox.getVerticalDeflectionFunction(fidelity=20), axisSameScale=False)
print(wingbox.getVerticalDeflectionAtY(AircraftProperties.Planform["span"] / 2))

wingbox.drawTwist(wingbox.getTwistFunction(), degrees=True)

from aircraftProperties import AircraftProperties
from Polygon import StringerType
from wingbox import Wingbox
import AerodynamicLoading
import numpy as np
import scipy as sp
from scipy import interpolate
from WingLoads import WingLoads
from WingLoads import PointForce



def getWingLoading(velocity, altitude, weight, loadFactor, engineThrustFactor, fuelMassDistribution=lambda y: 0, structuralMassDistribution=lambda y: 0, fuelMassFactor=1.0, neglectTangential=True):
    temp, p, rho = AerodynamicLoading.getISAParameters(altitude)

    S = AircraftProperties.Planform["surface area"]
    cL = (weight * loadFactor) / (1 / 2 * rho * (velocity**2) * S)
    forces = AerodynamicLoading.getNormalTangentialMomentAOA(cL, 0.375, velocity, altitude)

    aoa = forces[3]

    thrust = AircraftProperties.Engine["thrust"] * 1000 * engineThrustFactor

    # 7.369915250148878 is the chord at 10.02
    engine = PointForce([-7.369915250148878 * 0.375, 10.02, -1],
                        [-thrust, 0, -AircraftProperties.Engine["weight"]*loadFactor], "Engine")
    engineRotated = engine.getRotatedVectorAroundAxis(-aoa, [0, 1, 0])

    print(engineRotated.getMomentsAroundPoint([0, 10.02, 0]))

    #add fuel and structural weight
    normal = lambda y: forces[0](y) - (fuelMassDistribution(y) * fuelMassFactor + structuralMassDistribution(y)) * np.cos(aoa) * 9.80665 * loadFactor

    if neglectTangential:
        tangential = lambda y: 0
    else:
        tangential = lambda y: forces[1](y) - (fuelMassDistribution(y) * fuelMassFactor + structuralMassDistribution(y)) * np.sin(aoa) * 9.80665 * loadFactor

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

# amount stringers top, amount stringers bottom, endpoint
def wingboxLayoutHelper(sectionEndLocations, stringersTop, stringersBottom, stringerType, outerSparLocations, extraSpars, sparThicknesses):
    #section starting y values
    sectionStartingLocations = [0]
    sectionStartingLocations.extend(sectionEndLocations)

    #spar locations
    sparDistance = outerSparLocations[1] - outerSparLocations[0]

    #stringers
    stringerTypeTop = stringerType
    stringerTypeBottom = stringerType.getMirrorStringerX()
    topStringersOutput = []
    bottomStringersOutput = []
    #top stringers
    for sectionIndex in range(len(sectionEndLocations)):
        for stringerIndex in range(stringersTop[sectionIndex]):
            topStringersOutput.append([sectionStartingLocations[sectionIndex], sectionEndLocations[sectionIndex], outerSparLocations[0] + sparDistance / (stringersTop[sectionIndex] + 1) * (stringerIndex + 1), stringerTypeTop])
    #bottom stringers
    for sectionIndex in range(len(sectionEndLocations)):
        for stringerIndex in range(stringersBottom[sectionIndex]):
            bottomStringersOutput.append([sectionStartingLocations[sectionIndex], sectionEndLocations[sectionIndex], outerSparLocations[0] + sparDistance / (stringersBottom[sectionIndex] + 1) * (stringerIndex + 1), stringerTypeBottom])

    #spars
    sparThicknessFunc = sp.interpolate.interp1d(sectionEndLocations, sparThicknesses, kind='next', fill_value="extrapolate")
    sparsOutput = []
    sparsOutput.append([0, AircraftProperties.Planform["span"]/2+1, outerSparLocations[0], sparThicknessFunc])
    sparsOutput.append([0, AircraftProperties.Planform["span"]/2+1, outerSparLocations[1], sparThicknessFunc])
    for sectionIndex in range(len(sectionEndLocations)):
        for i in range(extraSpars[sectionIndex]):
            sparsOutput.append([sectionStartingLocations[sectionIndex], sectionEndLocations[sectionIndex], outerSparLocations[0] + (sparDistance/(extraSpars[sectionIndex] + 1)) * (i + 1), sparThicknessFunc])

    return topStringersOutput, bottomStringersOutput, sparsOutput, sparThicknessFunc
#wingbox definition start

scale = 200
stringerType = StringerType(
    [[0, 0], [0, -3 / scale], [17 / scale, -3 / scale], [17 / scale, -20 / scale], [20 / scale, -20 / scale],
     [20 / scale, 0]], [[0, 0], [20 / scale, 0]])

sections = [10, 15, 20, 30]

outerSparLocations = [0.15, 0.6]
extraSpars = [4, 0, 1, 0]
sparThicknesses = [0.1, 0.1, 0.1, 0.1]
flangeThicknesses = [0.1, 0.1, 0.1, 0.1]
stringersTop = [8, 6, 3, 2]
stringersBottom = [2, 2, 1, 0]

wingboxInputs = wingboxLayoutHelper(sectionEndLocations=sections, stringersTop=stringersTop, stringersBottom=stringersBottom, stringerType=stringerType, outerSparLocations=outerSparLocations, extraSpars=extraSpars, sparThicknesses=sparThicknesses)


yLocations = [-1, 30]
thicknesses = [-1, 0.01]
flangeThicknessTopFunc = sp.interpolate.interp1d(yLocations, thicknesses, kind='next', fill_value="extrapolate")

yLocations = [-1, 30]
thicknesses = [-1, 0.05]
flangeThicknessBottomFunc = sp.interpolate.interp1d(yLocations, thicknesses, kind='next', fill_value="extrapolate")

#wingbox Definition end



# rib location : fraction location of semispan, spars: [[absstart, absend, fractchordloc, thickness], []], sparThickness: [[[Thicknesses of spar], endpoint], [[Thicknesses of spar2], endpoint2]], stringersTop = [[starty, endy, fraction of chord, stringer]] flangeThickness = same as spar thickness
wingbox = Wingbox(ribLocations=sections, spars=wingboxInputs[2], stringersTop=wingboxInputs[0], stringersBottom=wingboxInputs[1], sparFlangeConnectionStringerShape=stringerType, flangeThicknesses=[flangeThicknessTopFunc, flangeThicknessBottomFunc], crosssectionAmount=200)
wingbox.draw(drawBottomStringers=False)


wingbox.drawStructuralMassDistribution()
wingbox.drawFuelMassDistribution()

wingloading = getWingLoading(velocity=232.412, altitude=AircraftProperties.Cruise_constants["cruise altitude"], weight=1866756, engineThrustFactor=0, loadFactor=1, fuelMassDistribution=wingbox.getFuelMassDistribution(), fuelMassFactor=0.8, neglectTangential=True)
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

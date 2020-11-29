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
    SF = AircraftProperties.Safety["ExternalLoadSF"]

    S = AircraftProperties.Planform["surface area"]
    cL = (weight * loadFactor * SF) / (1 / 2 * rho * (velocity**2) * S)
    forces = AerodynamicLoading.getNormalTangentialMomentAOA(cL, 0.375, velocity, altitude)
    aoa = 0

    thrust = AircraftProperties.Engine["thrust"] * 1000 * engineThrustFactor

    # 7.369915250148878 is the chord at 10.02
    engine = PointForce([-7.369915250148878 * 0.375, 10.02, -1], [-thrust, 0, -AircraftProperties.Engine["weight"]*loadFactor], "Engine")
    engineRotated = engine.getRotatedVectorAroundAxis(-aoa, [0, 1, 0])
    #add fuel and structural weight
    normal = lambda y: forces[0](y) - (fuelMassDistribution(y) * fuelMassFactor + structuralMassDistribution(y)) * np.cos(aoa) * 9.80665 * loadFactor

    if neglectTangential:
        engineRotated.forceVector[0] = 0
        tangential = lambda y: 0
    else:
        tangential = lambda y: forces[1](y) - (fuelMassDistribution(y) * fuelMassFactor + structuralMassDistribution(y)) * np.sin(aoa) * 9.80665 * loadFactor

    wingloading = WingLoads([engine], [tangential, None, normal], [None, forces[2], None])

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
def wingboxLayoutHelper(sectionEndLocations, stringersTop, stringersBottom, stringerType, outerSparLocations, extraSpars, sparThicknesses, flangeThicknessesTop, flangeThicknessesBottom):
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

    flangeThicknessTopFunc = sp.interpolate.interp1d(sectionEndLocations, flangeThicknessesTop, kind='next', fill_value="extrapolate")
    flangeThicknessBottomFunc = sp.interpolate.interp1d(sectionEndLocations, flangeThicknessesBottom, kind='next', fill_value="extrapolate")

    return topStringersOutput, bottomStringersOutput, sparsOutput, sparThicknessFunc, flangeThicknessTopFunc, flangeThicknessBottomFunc

def checkWingBox(loadingCases, wingbox):
    yieldStrength = AircraftProperties.WingboxMaterial["yield strength"]

    maxTwistDelfection = np.deg2rad(10)
    maxVerticalDeflection = AircraftProperties.Planform["span"]*0.15

    semispan = AircraftProperties.Planform["span"]/2

    fuelmassdistr = wingbox.getFuelMassDistribution()
    structmassdistr = wingbox.getStructuralMassDistribution()
    thrustlevels = [0.0, 1.0]

    for wingLoadingCase1 in loadingCases:
        for i in thrustlevels:
            loading = getWingLoading(velocity=wingLoadingCase1[0], altitude=wingLoadingCase1[2], weight=wingLoadingCase1[1], loadFactor=wingLoadingCase1[3], engineThrustFactor=i, fuelMassDistribution=fuelmassdistr, structuralMassDistribution=structmassdistr, fuelMassFactor=wingLoadingCase1[4], neglectTangential=True)
            wingbox.wingLoading = loading[0]

            # plot loading cases
            plotWingLoading(loading[0])

            #check normal stress
            max, min = wingbox.getMaximumMinimumNormalStress()
            if abs(max[0]) > yieldStrength or abs(min[0]) > yieldStrength:
                print("Max stress", max[0], "at y:", max[1].yLocation, "at point:", max[2])
                wingbox.drawCrosssection(max[1].yLocation, drawBendingStress=True)
                print("Min stress", min[0], "at y:", min[1].yLocation, "at point:", min[2])
                wingbox.drawCrosssection(min[1].yLocation, drawBendingStress=True)
                return False

            deflection = wingbox.getVerticalDeflectionAtY(semispan)
            print("Deflection:", deflection)
            if abs(deflection) > maxVerticalDeflection:
                print("Exceeded Vertical Deflection Limit of", maxVerticalDeflection, "m. Deflection is:", deflection, "m.")
                wingbox.drawDeflection(wingbox.getVerticalDeflectionFunction())
                return False

            twist = wingbox.getTwist(semispan)
            print("Twist:", np.rad2deg(twist))
            if abs(twist) > maxTwistDelfection:
                print("Exceeded Twist Deflection Limit of", np.rad2deg(maxTwistDelfection), "deg. Deflection is:", np.rad2deg(twist), "deg.")
                wingbox.drawTwist(wingbox.getTwistFunction(), degrees=True)
                return False
    return True

v =114.60056632591561
weight = 1539078
altitude = 0
loadFactor1 = 2.5
fuelFactor = 0.0

loadFactor2 = -1

thrustlevels = [0.0, 1.0]

loadCases = [[v, weight, altitude, loadFactor1, fuelFactor], [v, weight, altitude, loadFactor2, fuelFactor]]

#wingbox definition start

scale = 4 # variable

stringerSize = 1000 / scale # 1000, 20x20mm, 2mm thick   * scale
stringerType = StringerType(
    [[0, 0], [0, -2 / stringerSize], [18 / stringerSize, -2 / stringerSize], [18 / stringerSize, -20 / stringerSize], [20 / stringerSize, -20 / stringerSize],
     [20 / stringerSize, 0]], [[0, 0], [20 / stringerSize, 0]])


outerSparLocations = [0.15, 0.6] # variable

sectionEnds =               [2.5,   5,      8,      11,     15,     19.5,     24,   30] # m
#per section
extraSpars =                [0,     0,      0,      0,      0,      0,      0,      0] # m
sparThicknesses =           [0.014, 0.013,  0.011,  0.009,  0.010,  0.008,   0.007,  0.005]  # m
flangeThicknessesTop =      [0.016, 0.015,  0.014,  0.013,  0.011,  0.009,  0.008,  0.006] # m
flangeThicknessesBottom =   [0.016, 0.015,  0.014,  0.013,  0.011,  0.009,  0.008,  0.006] # m
stringersTop =              [34,    32,     30,     17,     14,     5,      2,      1] # m
stringersBottom =           [25,    22,     20,     13,     9,      3,      2,      1] # m

wingboxInputs = wingboxLayoutHelper(sectionEndLocations=sectionEnds, stringersTop=stringersTop, stringersBottom=stringersBottom, stringerType=stringerType, outerSparLocations=outerSparLocations, extraSpars=extraSpars, sparThicknesses=sparThicknesses, flangeThicknessesTop=flangeThicknessesTop, flangeThicknessesBottom=flangeThicknessesBottom)


wingbox = Wingbox(ribLocations=sectionEnds, spars=wingboxInputs[2], stringersTop=wingboxInputs[0], stringersBottom=wingboxInputs[1], sparFlangeConnectionStringerShape=stringerType, flangeThicknesses=[wingboxInputs[4], wingboxInputs[5]], crosssectionAmount=200)


# draw wingbox
wingbox.draw(drawBottomStringers=False)
wingbox.drawCrosssection(2.5, drawCentroid=True, drawSidewallCenterlines=True)
#wingbox.drawInertias()

if checkWingBox(loadCases, wingbox):
    print("YAY it worked")
    print("Mass:", wingbox.totalMass, "kg")
else:
    print("Doesnt work")
    print("Mass:", wingbox.totalMass, "kg")


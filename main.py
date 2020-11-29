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

    #aoa = forces[3]
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

            #plot loading cases
            #plotWingLoading(loading[0])

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
                print("Exceeded Twist Deflection Limit of", np.rad2deg(maxVerticalDeflection), "deg. Deflection is:", np.rad2deg(deflection), "deg.")
                wingbox.drawTwist(wingbox.getTwistFunction())
                return False
    return True

v =114.60056632591561
weight = 1539078
altitude = 0
loadFactor1 = 2.5
fuelFactor = 0.0

loadFactor2 = -1.0

loadCases = [[v, weight, altitude, loadFactor1, fuelFactor], [v, weight, altitude, loadFactor2, fuelFactor]]

#wingbox definition start

scale = 2 # variable

stringerSize = 1000 / scale # 1000, 20x20mm, 1.5mm thick
stringerType = StringerType(
    [[0, 0], [0, -1.5 / stringerSize], [18.5 / stringerSize, -1.5 / stringerSize], [18.5 / stringerSize, -20 / stringerSize], [20 / stringerSize, -20 / stringerSize],
     [20 / stringerSize, 0]], [[0, 0], [20 / stringerSize, 0]])

sectionEnds = [5, 10.02, 15, 20, 30] # variable

outerSparLocations = [0.15, 0.6] # variable

#per section
extraSpars = [1, 1, 1, 0, 0] # variable
sparThicknesses = [0.01, 0.001, 0.01, 0.01, 0.01]  # variable
flangeThicknessesTop = [0.01, 0.01, 0.01, 0.01, 0.01] # variable
flangeThicknessesBottom = [0.01, 0.01, 0.01, 0.01, 0.01] # variable
stringersTop = [35, 30, 25, 15, 10] # variable
stringersBottom = [25, 15, 10, 5, 5] # variable

wingboxInputs = wingboxLayoutHelper(sectionEndLocations=sectionEnds, stringersTop=stringersTop, stringersBottom=stringersBottom, stringerType=stringerType, outerSparLocations=outerSparLocations, extraSpars=extraSpars, sparThicknesses=sparThicknesses, flangeThicknessesTop=flangeThicknessesTop, flangeThicknessesBottom=flangeThicknessesBottom)

wingbox = Wingbox(ribLocations=sectionEnds, spars=wingboxInputs[2], stringersTop=wingboxInputs[0], stringersBottom=wingboxInputs[1], sparFlangeConnectionStringerShape=stringerType, flangeThicknesses=[wingboxInputs[4], wingboxInputs[5]], crosssectionAmount=200)

#draw wingbox
#wingbox.draw(drawBottomStringers=False)
wingbox.drawCrosssection(2.5, drawCentroid=True, drawSidewallCenterlines=True)
#wingbox.drawInertias()

if checkWingBox(loadCases, wingbox):
    print("YAY it worked")
    print("Mass:", wingbox.totalMass, "kg")
else:
    print("Doesnt work")


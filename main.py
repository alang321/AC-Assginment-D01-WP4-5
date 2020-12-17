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

    if loadFactor > 0:
        cL = (weight * loadFactor * SF) / (1 / 2 * rho * (velocity**2) * S)
        forces = AerodynamicLoading.getNormalTangentialMomentAOA(cL, 0.375, velocity, altitude)
    else:
        cL = (weight * loadFactor * -1 * SF) / (1 / 2 * rho * (velocity**2) * S)
        invertedForces = AerodynamicLoading.getNormalTangentialMomentAOA(cL, 0.375, velocity, altitude)
        forces = []
        for i in range(3):
            forces.append(lambda y: invertedForces[i](y) * -1)
            forces.append(invertedForces[3] * -1)
    aoa = 0

    thrust = AircraftProperties.Engine["thrust"] * 1000 * engineThrustFactor

    # 7.369915250148878 is the chord at 10.02
    engine = PointForce([-7.369915250148878 * 0.375, 10.02, -1], [-thrust, 0, -AircraftProperties.Engine["weight"]*loadFactor], "Engine")
    engineRotated = engine.getRotatedVectorAroundAxis(-aoa, [0, 1, 0])
    #add fuel and structural weight
    normal = lambda y: forces[0](y) - (fuelMassDistribution(y) * fuelMassFactor + structuralMassDistribution(y)) * np.cos(aoa) * 9.80665 * loadFactor

    if neglectTangential:
        tangential = lambda y: 0
    else:
        tangential = lambda y: forces[1](y) - (fuelMassDistribution(y) * fuelMassFactor + structuralMassDistribution(y)) * np.sin(aoa) * 9.80665 * loadFactor

    wingloading = WingLoads([engine], [tangential, None, normal], [None, forces[2], None])

    return wingloading, aoa

def plotWingLoading(loading):
    print("External Forces")
    loading.drawExternalForces(2)

    print("External Moments")
    loading.drawExternalMoments(1)

    print("Shear Forces")
    loading.drawShearForce(2)

    print("Internal Moments")
    loading.drawInternalMoment(0)
    loading.drawInternalMoment(1)

# amount stringers top, amount stringers bottom, endpoint
def wingboxLayoutHelper(sectionEndLocations, stringersTop, stringersBottom, stringerType, sparCapSide, sparCapCenter, outerSparLocations, extraSpars, sparThicknesses, flangeThicknessesTop, flangeThicknessesBottom, ribThickness):
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
        cells = extraSpars[sectionIndex] + 1

        extra = stringersTop[sectionIndex] % cells

        stringersPerCell = []
        for i in range(cells):
            stringersPerCell.append(stringersTop[sectionIndex]//cells)

            if i < extra:
                stringersPerCell[-1] += 1

        for i in range(cells):
            chordAtMedian = 9.77 * (1 + ((2 * (0.2996 - 1)) / 57.1362641792155) * (sectionStartingLocations[sectionIndex] + sectionEndLocations[sectionIndex])/2)

            frontlength = sparCapCenter.baseLength
            aftlength = sparCapCenter.baseLength
            if i == 0:
                frontlength = sparCapSide.baseLength
            if i == cells - 1:
                aftlength = sparCapSide.baseLength

            start = outerSparLocations[0] + (sparDistance/(extraSpars[sectionIndex] + 1)) * i + frontlength/chordAtMedian - stringerTypeTop.baseLength/2/chordAtMedian
            end = outerSparLocations[0] + (sparDistance/(extraSpars[sectionIndex] + 1)) * (i + 1) - aftlength/chordAtMedian + stringerTypeTop.baseLength/2/chordAtMedian

            for stringerIndex in range(stringersPerCell[i]):
                topStringersOutput.append([sectionStartingLocations[sectionIndex], sectionEndLocations[sectionIndex], start + (end - start)/(stringersPerCell[i] + 1) * (stringerIndex + 1), stringerTypeTop])
    #bottom stringers
    for sectionIndex in range(len(sectionEndLocations)):
        cells = extraSpars[sectionIndex] + 1

        extra = stringersBottom[sectionIndex] % cells

        stringersPerCell = []
        for i in range(cells):
            stringersPerCell.append(stringersBottom[sectionIndex] // cells)

            if i < extra:
                stringersPerCell[-1] += 1

        for i in range(cells):
            chordAtMedian = 9.77 * (1 + ((2 * (0.2996 - 1)) / 57.1362641792155) * (
                        sectionStartingLocations[sectionIndex] + sectionEndLocations[sectionIndex]) / 2)

            frontlength = sparCapCenter.baseLength
            aftlength = sparCapCenter.baseLength
            if i == 0:
                frontlength = sparCapSide.baseLength
            if i == cells - 1:
                aftlength = sparCapSide.baseLength

            start = outerSparLocations[0] + (sparDistance / (extraSpars[sectionIndex] + 1)) * i + frontlength / chordAtMedian - stringerTypeBottom.baseLength/2/chordAtMedian
            end = outerSparLocations[0] + (sparDistance / (extraSpars[sectionIndex] + 1)) * (i + 1) - aftlength / chordAtMedian + stringerTypeBottom.baseLength/2/chordAtMedian

            for stringerIndex in range(stringersPerCell[i]):
                bottomStringersOutput.append([sectionStartingLocations[sectionIndex], sectionEndLocations[sectionIndex],start + (end - start) / (stringersPerCell[i] + 1) * (stringerIndex + 1), stringerTypeBottom])
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

    return Wingbox(ribLocations=[0, *sectionEndLocations], spars=sparsOutput, stringersTop=topStringersOutput, stringersBottom=bottomStringersOutput, sparCapSide=sparCapSide, sparCapCenter=sparCapCenter, ribThickness=ribThickness, flangeThicknesses=[flangeThicknessTopFunc, flangeThicknessBottomFunc], crosssectionAmount=200)

def checkWingBox(loadingCases, wingbox):
    yieldStrength = AircraftProperties.WingboxMaterial["yield strength"]
    edgecrackStrength = AircraftProperties.WingboxMaterial["edge crack strength"]
    centercrackStrength = AircraftProperties.WingboxMaterial["center crack strength"]

    maxTwistDelfection = np.deg2rad(10)
    maxVerticalDeflection = AircraftProperties.Planform["span"]*0.15

    semispan = AircraftProperties.Planform["span"]/2

    fuelmassdistr = wingbox.getFuelMassDistribution()
    structmassdistr = wingbox.getStructuralMassDistribution()
    thrustlevels = [1.0]

    for wingLoadingCase1 in loadingCases:
        for i in thrustlevels:
            loading = getWingLoading(velocity=wingLoadingCase1[0], altitude=wingLoadingCase1[2], weight=wingLoadingCase1[1], loadFactor=wingLoadingCase1[3], engineThrustFactor=i, fuelMassDistribution=fuelmassdistr, structuralMassDistribution=structmassdistr, fuelMassFactor=wingLoadingCase1[4], neglectTangential=True)
            wingbox.wingLoading = loading[0]

            # plot loading cases
            #plotWingLoading(loading[0])

            wingbox.drawNormalStressSafetyMargin()

            wingbox.drawFlangeBucklingSafetyMargin()
            wingbox.drawShearWebBucklingSafetyMargin()
            wingbox.drawColumnBucklingSafetyMargin()
            wingbox.drawInterRivetBucklingSafetyMargin()

            print(wingbox.jFuncY(28.06))

            if wingbox.checkFlangeBuckling():
                return False

            if wingbox.checkShearWebBuckling():
                return False

            if wingbox.checkColumnBuckling():
                return False

            if wingbox.checkInterRivetBuckling():
                return False

            if wingbox.jFuncY(28.06) <= 7.7 * 10**-4:
                print("Aileron reversal occurs. J is", str(wingbox.jFuncY(28.06)), "at", str(28.06), "[m] but should be atleast", str(7.7 * 10**-4))
                return False

            #check normal stress
            max, min = wingbox.getMaximumMinimumNormalStress()
            if abs(max[0]) > yieldStrength or abs(min[0]) > yieldStrength:
                print("Max stress", max[0], "at y:", max[1].yLocation, "at point:", max[2])
                wingbox.drawCrosssection(max[1].yLocation, drawBendingStress=True)
                print("Min stress", min[0], "at y:", min[1].yLocation, "at point:", min[2])
                wingbox.drawCrosssection(min[1].yLocation, drawBendingStress=True)
                return False

            #check with cracks
            if abs(max[0]) > edgecrackStrength:
                print("Edge cracks will lead to failure.")
                return False
            if abs(max[0]) > centercrackStrength:
                print("Center cracks will lead to failure.")
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
fuelFactor = 0.0

loadFactor1 = 2.5
loadFactor2 = -1

loadCases = [[v, weight, altitude, loadFactor1, fuelFactor], [v, weight, altitude, loadFactor2, fuelFactor]]

#wingbox definition start


scale = 8
sideCap = StringerType([[0, 0], [0, -3], [20, -3], [20, -23], [23, -23], [23, 0]], [[0, 0], [23, 0]], rivetPoints=[[3, 0]], isSparCap=True)
sideCap = sideCap.getScaledStringer(1/1000 * scale)

scale = 8
centerCap = StringerType([[0, 0], [0, -1.5], [21.5, -1.5], [21.5, -23], [23, -23], [23, 0]], [[0, 0], [23, 0]], rivetPoints=[[20, 0]], isSparCap=True)
centerCap = centerCap.getScaledStringer(1/1000 * scale)

scale = 3.7
extrudedt = StringerType([[0, 0], [0, -1.5], [8.5, -1.5], [8.5, -18.5], [3, -18.5], [3, -20], [15.5, -20], [15.5, -18.5], [10, -18.5], [10, -1.5], [18.5, -1.5], [18.5, 0]],       [[0, 0], [18.5, 0]])
extrudedt = extrudedt.getScaledStringer((1/1000)*scale)

outerSparLocations = [0.15, 0.6] # variable

sectionEnds =               [1.05,    2.15,    3.35,   4.6,    5.9,      7.25,     8.6,      10,       11.4,    12.8,   14.2,   15.6,   17,     18.4,   19.9,   21.45,  23,     24.5,   26.5,   28.58] # m
#per section
extraSpars =                [1,       1,       1,      1,      1,        1,        1,        1,        0,       0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0] # m
sparThicknesses =           [0.0092,  0.0091,  0.009,  0.0089, 0.0086,   0.0084,   0.0081,   0.0078,   0.0082,  0.0078, 0.0074, 0.0068, 0.0066, 0.006,  0.0054, 0.0048, 0.0042, 0.0038, 0.0032, 0.0028]  # m
flangeThicknessesTop =      [0.018,   0.0163,  0.015,  0.0141, 0.0131,   0.0122,   0.0111,   0.0105,   0.0115,  0.010,  0.0096, 0.0087, 0.0078, 0.0069, 0.0062, 0.0058, 0.0053, 0.0046, 0.0043, 0.0038] # m
flangeThicknessesBottom =   [0.017,   0.0153,  0.014,  0.0131, 0.0121,   0.0112,   0.0101,   0.0095,   0.0105,  0.009,  0.0086, 0.0072, 0.0067, 0.0059, 0.0045, 0.004,  0.0035, 0.0033, 0.0029, 0.0028] # m
stringersTop =              [20,      18,      17,     15,     14,       13,       12,       12,       10,      10,     9,      9,      8,      8,      7,      7,      6,      5,      4,      2] # m
stringersBottom =           [16,      15,      14,     13,     12,       10,       10,       10,       8,       6,      5,      4,      4,      3,      3,      3,      2,      2,      2,      1] # m

ribThickness = 0.003 # m, set this to the lowest skin thickness value


wingbox = wingboxLayoutHelper(sectionEndLocations=sectionEnds, stringersTop=stringersTop, stringersBottom=stringersBottom, stringerType=extrudedt, sparCapSide=sideCap, sparCapCenter=centerCap, outerSparLocations=outerSparLocations, extraSpars=extraSpars, sparThicknesses=sparThicknesses, flangeThicknessesTop=flangeThicknessesTop, flangeThicknessesBottom=flangeThicknessesBottom, ribThickness=ribThickness)

wingbox.drawMinimumRivetPitch()
wingbox.drawFlangeThickness()
wingbox.drawCombinedSparThickness()

print(wingbox.totalMass)
print("Fuel Volume:", wingbox.internalVolume)

# draw wingbox
wingbox.draw(drawBottomStringers=False)
#wingbox.drawCrosssection(2.5, drawCentroid=True, drawSidewallCenterlines=True)
#wingbox.getGeneratedCrosssectionAtY(2.5).getMaxShearStress(0, 0)
#wingbox.drawInertias()

if checkWingBox(loadCases, wingbox):
    print("YAY it worked")
    print("Mass:", wingbox.totalMass, "kg")
    print("Fuel Volume:", wingbox.internalVolume)
else:
    print("Doesnt work")
    print("Mass:", wingbox.totalMass, "kg")
    print("Fuel Volume:", wingbox.internalVolume)


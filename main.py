from aircraftProperties import AircraftProperties
from Polygon import StringerType
from wingbox import Wingbox
import AerodynamicLoading
import numpy as np
import scipy as sp
from scipy import interpolate
from WingLoads import WingLoads
from WingLoads import PointForce

AerodynamicLoading.drawAerodynamicCoefficients()

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

    for wingLoadingCase in loadingCases:

        wingbox.wingLoading = wingLoadingCase

        #check normal stress
        max, min = wingbox.getMaximumorMinimumNormalStress()
        if abs(max[0]) > yieldStrength or abs(min[0]) > yieldStrength:
            print("Max stress", max[0], "at y:", max[1].yLocation, "at point:", max[2])
            wingbox.drawCrosssection(max[1].yLocation, drawBendingStress=True)
            print("Min stress", min[0], "at y:", min[1].yLocation, "at point:", min[2])
            wingbox.drawCrosssection(min[1].yLocation, drawBendingStress=True)
            return False

        deflection = wingbox.getVerticalDeflectionAtY(semispan)
        if abs(deflection) > maxVerticalDeflection:
            print("Exceeded Vertical Deflection Limit of", maxVerticalDeflection, "m. Deflection is:", deflection, "m.")
            wingbox.drawDeflection(wingbox.getVerticalDeflectionFunction())
            return False

        twist = wingbox.getTwist(semispan)
        if abs(twist) > maxTwistDelfection:
            print("Exceeded Twist Deflection Limit of", np.rad2deg(maxVerticalDeflection), "deg. Deflection is:", np.rad2deg(deflection), "deg.")
            wingbox.drawTwist(wingbox.getTwistFunction())
            return False
    return True

#fuel mass distribution
fuelDistrList = [904.2902174010153, 897.8569544055961, 891.4466560950258, 885.0593224693056, 878.6949535284349, 872.3535492724135, 866.0351097012415, 859.7396348149188, 853.4671246134459, 847.2175790968226, 840.9909982650483, 834.7873821181239, 828.6067306560492, 822.4490438788235, 816.3143217864475, 810.2025643789209, 804.1137716562439, 798.0479436184164, 792.0050802654383, 785.9851815973099, 779.9882476140308, 774.0142783156007, 768.0632737020209, 762.1352337732902, 756.2301585294093, 750.3480479703773, 744.4889020961954, 738.6527209068626, 732.8395044023791, 727.0492525827456, 721.2819654479614, 715.5376429980266, 709.8162852329413, 704.1178921527053, 698.4424637573193, 692.7900000467823, 687.1605010210949, 681.5539666802572, 675.9703970242688, 670.4097920531299, 664.8721517668405, 659.3574761654006, 653.8657652488105, 648.3970190170694, 642.9512374701778, 637.5284206081361, 632.1285684309435, 626.7516809386004, 621.397758131107, 616.0668000084631, 610.7588065706683, 605.4737778177234, 600.2117137496278, 594.9726143663817, 589.7564796679851, 584.563309654438, 579.3931043257406, 574.2458636818924, 569.1215877228935, 564.0202764487444, 558.9419298594448, 553.8865479549946, 548.854130735394, 543.8446782006429, 538.8581903507411, 533.8946671856888, 528.954108705486, 524.0365149101328, 519.141885799629, 514.2702213739748, 509.42152163316996, 504.5957865772148, 499.79301620610875, 495.0132105198525, 490.2563695184456, 485.5224932018882, 480.8115815701803, 476.12363462332183, 471.45865236131306, 466.81663478415356, 462.1975818918435, 457.6014936843831, 453.0283701617722, 448.47821132401066, 443.9510171710987, 439.4467877030363, 434.9655229198232, 430.50722282145966, 426.07188740794567, 421.65951667928107, 417.270110635466, 412.9036692765003, 408.56019260238406, 404.2396806131176, 399.9421333087004, 395.66755068913284, 391.41593275441454, 387.187279504546, 382.98159093952677, 378.798867059357, 374.63910786403693, 370.5023133535662, 366.38848352794486, 362.2976183871732, 358.2297179312508, 355.47660738055987, 351.44721584241273, 347.440788989115, 343.4573268206667, 339.49682933706794, 335.5592965383187, 331.6447284244189, 327.7531249953686, 323.8844862511677, 320.0388121918164, 316.21610281731455, 312.41635812766214, 308.6395781228593, 304.8857628029058, 301.1549121678019, 297.4470262175476, 293.76210495214275, 290.1001483715871, 286.4611564758811, 282.8451292650246, 279.2520667390178, 275.68196889786026, 272.1348357415521, 268.6106672700937, 265.1094634834846, 261.63122438172496, 258.175949964815, 254.74364023275427, 251.3342951855434, 247.9479148231816, 244.58449914566933, 241.24404815300684, 237.92656184519365, 234.63204022222996, 231.36048328411576, 228.11189103085104, 224.88626346243606, 221.68360057887028, 218.50390238015387, 215.34716886628712, 212.21340003726985, 209.10259589310218, 206.01475643378373, 202.949881659315, 199.90797156969566, 196.8890261649259, 193.8930454450054, 190.92002940993459, 187.9699780597133, 185.04289139434132, 182.1387694138189, 179.257612118146, 176.3994195073226, 173.56419158134855, 170.75192834022408, 167.9626297839492, 165.1962959125237, 162.4529267259476, 159.73252222422116, 157.03508240734413, 154.3606072753166, 151.70909682813843, 149.08055106580989, 146.47496998833077, 143.89235359570125, 141.33270188792102, 138.7960148649904, 136.28229252690932, 133.79153487367768, 131.32374190529546, 128.87891362176276, 126.45705002307959, 124.05815110924588, 121.68221688026162, 119.32924733612681, 116.99924247684163, 114.69220230240589, 112.40812681281955, 110.14701600808272, 107.90886988819544, 105.6936884531577, 103.5014717029693, 101.33221963763037, 99.18593225714099, 97.06260956150122, 94.96225155071079, 92.88485822476983, 90.83042958367841, 88.79896562743647, 86.79046635604404, 84.80493176950102, 82.8423618678076, 80.9027566509636, 78.98611611896911, 77.09244027182403]
FuelYLocations = [0.0, 0.14355845271159673, 0.28711690542319346, 0.4306753581347902, 0.5742338108463869, 0.7177922635579836, 0.8613507162695804, 1.004909168981177, 1.1484676216927738, 1.2920260744043706, 1.4355845271159673, 1.579142979827564, 1.7227014325391607, 1.8662598852507575, 2.009818337962354, 2.153376790673951, 2.2969352433855477, 2.4404936960971444, 2.584052148808741, 2.727610601520338, 2.8711690542319346, 3.0147275069435313, 3.158285959655128, 3.3018444123667248, 3.4454028650783215, 3.5889613177899182, 3.732519770501515, 3.8760782232131117, 4.019636675924708, 4.163195128636305, 4.306753581347902, 4.450312034059499, 4.593870486771095, 4.737428939482692, 4.880987392194289, 5.0245458449058855, 5.168104297617482, 5.311662750329079, 5.455221203040676, 5.598779655752272, 5.742338108463869, 5.885896561175466, 6.029455013887063, 6.173013466598659, 6.316571919310256, 6.460130372021853, 6.6036888247334495, 6.747247277445046, 6.890805730156643, 7.03436418286824, 7.1779226355798365, 7.321481088291433, 7.46503954100303, 7.608597993714627, 7.752156446426223, 7.89571489913782, 8.039273351849417, 8.182831804561014, 8.32639025727261, 8.469948709984207, 8.613507162695804, 8.7570656154074, 8.900624068118997, 9.044182520830594, 9.18774097354219, 9.331299426253787, 9.474857878965384, 9.61841633167698, 9.761974784388578, 9.905533237100174, 10.049091689811771, 10.192650142523368, 10.336208595234964, 10.479767047946561, 10.623325500658158, 10.766883953369755, 10.910442406081351, 11.054000858792948, 11.197559311504545, 11.341117764216142, 11.484676216927738, 11.628234669639335, 11.771793122350932, 11.915351575062529, 12.058910027774125, 12.202468480485722, 12.346026933197319, 12.489585385908915, 12.633143838620512, 12.776702291332109, 12.920260744043706, 13.063819196755302, 13.207377649466899, 13.350936102178496, 13.494494554890093, 13.63805300760169, 13.781611460313286, 13.925169913024883, 14.06872836573648, 14.212286818448076, 14.355845271159673, 14.49940372387127, 14.642962176582866, 14.786520629294463, 14.93007908200606, 15.073637534717657, 15.217195987429253, 15.36075444014085, 15.504312892852447, 15.647871345564043, 15.79142979827564, 15.934988250987237, 16.078546703698834, 16.22210515641043, 16.365663609122027, 16.509222061833626, 16.65278051454522, 16.796338967256816, 16.939897419968414, 17.083455872680013, 17.227014325391607, 17.370572778103202, 17.5141312308148, 17.6576896835264, 17.801248136237994, 17.94480658894959, 18.088365041661188, 18.231923494372786, 18.37548194708438, 18.519040399795976, 18.662598852507575, 18.806157305219173, 18.94971575793077, 19.093274210642363, 19.23683266335396, 19.38039111606556, 19.523949568777155, 19.66750802148875, 19.81106647420035, 19.954624926911947, 20.098183379623542, 20.241741832335137, 20.385300285046736, 20.528858737758334, 20.67241719046993, 20.815975643181524, 20.959534095893122, 21.10309254860472, 21.246651001316316, 21.39020945402791, 21.53376790673951, 21.677326359451108, 21.820884812162703, 21.964443264874298, 22.108001717585896, 22.251560170297495, 22.39511862300909, 22.538677075720685, 22.682235528432283, 22.82579398114388, 22.969352433855477, 23.11291088656707, 23.25646933927867, 23.40002779199027, 23.543586244701864, 23.68714469741346, 23.830703150125057, 23.974261602836656, 24.11782005554825, 24.261378508259845, 24.404936960971444, 24.548495413683042, 24.692053866394637, 24.835612319106232, 24.97917077181783, 25.12272922452943, 25.266287677241024, 25.40984612995262, 25.553404582664218, 25.696963035375816, 25.84052148808741, 25.984079940799006, 26.127638393510605, 26.271196846222203, 26.414755298933798, 26.558313751645393, 26.70187220435699, 26.84543065706859, 26.988989109780185, 27.13254756249178, 27.27610601520338, 27.419664467914977, 27.563222920626572, 27.706781373338167, 27.850339826049765, 27.993898278761364, 28.13745673147296, 28.281015184184554, 28.424573636896152, 28.56813208960775, ]
fuelMassFunction = sp.interpolate.interp1d(FuelYLocations, fuelDistrList)

wingLoadingMaximum = 0
wingLoadingMinimum = 0

wingloading = getWingLoading(velocity=232.412, altitude=AircraftProperties.Cruise_constants["cruise altitude"], weight=1866756, engineThrustFactor=0, loadFactor=1, fuelMassDistribution=wingbox.getFuelMassDistribution(), fuelMassFactor=0.8, neglectTangential=True)


#wingbox definition start

scale = 200
stringerType = StringerType(
    [[0, 0], [0, -1.5 / scale], [18.5 / scale, -1.5 / scale], [18.5 / scale, -20 / scale], [20 / scale, -20 / scale],
     [20 / scale, 0]], [[0, 0], [20 / scale, 0]])

sectionEnds = [10, 15, 20, 30]

outerSparLocations = [0.15, 0.6]
extraSpars = [1, 1, 0, 0]
sparThicknesses = [0.01, 0.01, 0.01, 0.01]
flangeThicknessesTop = [0.01, 0.01, 0.01, 0.01]
flangeThicknessesBottom = [0.01, 0.01, 0.01, 0.01]
stringersTop = [6, 4, 3, 2]
stringersBottom = [2, 2, 1, 0]

wingboxInputs = wingboxLayoutHelper(sectionEndLocations=sectionEnds, stringersTop=stringersTop, stringersBottom=stringersBottom, stringerType=stringerType, outerSparLocations=outerSparLocations, extraSpars=extraSpars, sparThicknesses=sparThicknesses, flangeThicknessesTop=flangeThicknessesTop, flangeThicknessesBottom=flangeThicknessesBottom)

wingbox = Wingbox(ribLocations=sectionEnds, spars=wingboxInputs[2], stringersTop=wingboxInputs[0], stringersBottom=wingboxInputs[1], sparFlangeConnectionStringerShape=stringerType, flangeThicknesses=[wingboxInputs[4], wingboxInputs[5]], crosssectionAmount=200)

wingbox.draw(drawBottomStringers=False)
wingbox.drawCrosssection(5, drawCentroid=True, drawSidewallCenterlines=True)
wingbox.drawInertias()


#wingbox Definition end



# rib location : fraction location of semispan, spars: [[absstart, absend, fractchordloc, thickness], []], sparThickness: [[[Thicknesses of spar], endpoint], [[Thicknesses of spar2], endpoint2]], stringersTop = [[starty, endy, fraction of chord, stringer]] flangeThickness = same as spar thickness



#wingbox.drawStructuralMassDistribution()
#wingbox.drawFuelMassDistribution()

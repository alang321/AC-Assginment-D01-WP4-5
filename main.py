from aircraftProperties import AircraftProperties
from wingboxCrosssection import WingboxCrossection
from Polygon import StringerType
from wingbox import Wingbox
import AerodynamicLoading
import numpy as np
from distributedLoad import distributedLoad



#aerodynamicLoadingDistribution.drawAerodynamicCoefficients()

v = 232.412
h = AircraftProperties.Cruise_constants["cruise altitude"]
rho = AircraftProperties.Cruise_constants["density at cruise"]
S = 362.73
cL = 1866756/((1/2) * rho * v**2 * S)
forces = AerodynamicLoading.getNormalTangentialMomentAOA(cL, 0.375, v, h)
print(np.rad2deg(forces[3]))

AerodynamicLoading.drawAerodynamicCoefficients()
AerodynamicLoading.drawAerodynamicForces(forces)


#Add normal force
normalForce = distributedLoad([0, AircraftProperties.Planform["span"] / 2])

normalForce.addDistibutedLoad(forces[0])
normalForce.addPointLoad(-AircraftProperties.Engine["weight"], 10, "Engine")

normalForce.drawLoads(fidelity=100)
normalForce.getShearForce(fidelity=60)
normalForce.getMomentDistribution(fidelity=20)
normalForce.drawShear(fidelity=100)

#tangentialForce = distributedLoad([0, AircraftProperties.Planform["span"] / 2])
#tangentialForce.addDistibutedLoad(forces[1])
#tangentialForce.getShearForce(fidelity=20)
#tangentialForce.getMomentDistribution(fidelity=20)

twistMoment = distributedLoad([0, AircraftProperties.Planform["span"] / 2])
twistMoment.addMomentDistribution(forces[2])
twistMoment.drawMoment()




stringer = StringerType([[0, 0], [0, -3/500], [17/500, -3/200], [17/500, -20/500], [20/500, -20/500], [20/500, 0]], [[0, 0], [20/500, 0]])
stringerTop = stringer
stringerBottom = stringer.getMirrorStringerX()


# rib location : fraction location of semispan, spars: [[absstart, absend, fractchordloc, thickness], []], sparThickness: [[[Thicknesses of spar], endpoint], [[Thicknesses of spar2], endpoint2]], stringersTop = [[starty, endy, fraction of chord, stringer]] flangeThickness = same as spar thickness
wingbox = Wingbox(ribLocations=[], spars=[[0, 30, 0.15, 0.01], [0, 30, 0.6, 0.01], [0, 15, 0.375, 0.01]], stringersTop=[[0, 35, 0.3, stringer]], stringersBottom=[[0, 10, 0.4, stringerBottom]], sparFlangeConnectionStringerShape=stringer, flangeThicknesses=[[[0.01, 0.01], 30]], crosssectionAmount=200)
#wingbox = Wingbox(ribLocations=[], spars=[[0, 30, 0.16, 0.01], [0, 30, 0.65, 0.01], [0, 15, 0.405, 0.01]], stringersTop=[[0, 35, 0.3, stringer]], stringersBottom=[[0, 10, 0.4, stringerBottom]], sparFlangeConnectionStringerShape=stringer, flangeThicknesses=[[[0.01, 0.01], 30]], crosssectionAmount=200)
wingbox.draw(top=False)

#set moments
wingbox.shearzx = normalForce.getShearForce(fidelity=20)
wingbox.mxFunc = normalForce.getMomentDistribution()
wingbox.myFunc = twistMoment.getMomentDistribution()

wingbox.drawInertias()
wingbox.drawMaterialMassDistribution()
wingbox.drawFuelMassDistribution()

print(wingbox.getGeneratedCrosssectionAtY(0).getMaxShearStress(50000, 100000))

posY = 5
wingbox.getGeneratedCrosssectionAtY(posY).drawCrosssection(drawCentroid=True, drawSidewallCenterlines=True, drawBendingStress=True, Mx=wingbox.mxFunc(posY), Mz=0)

max, min = wingbox.getMaximumMinimumNormalStress()
print("Max stress", max[0], "at y:", max[1].yLocation, "at point:", max[2])
max[1].drawCrosssection(drawBendingStress=True, Mx=wingbox.mxFunc(max[1].yLocation), Mz=0)
print("Min stress", min[0], "at y:", min[1].yLocation, "at point:", min[2])
min[1].drawCrosssection(drawBendingStress=True, Mx=wingbox.mxFunc(min[1].yLocation), Mz=0)


#deflection
wingbox.drawDeflection(wingbox.getVerticalDeflectionFunction(fidelity=20), axisSameScale=False)
print(wingbox.getVerticalDeflectionAtY(AircraftProperties.Planform["span"] / 2))

wingbox.drawTwist(wingbox.getTwistFunction(), degrees=True)

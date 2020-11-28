import matplotlib.pyplot as plt
from aircraftProperties import AircraftProperties
from Polygon import Polygon
from wingboxCrosssection import WingboxCrossection
from scipy.interpolate import interp1d
from scipy.integrate import quad
from WingLoads import WingLoads
import numpy as np

class Wingbox:
    __fuelVolumeReq = AircraftProperties.Fuel["total fuel volume required"] / 2
    __fuelDensity = AircraftProperties.Fuel["fuel density"]
    __density = AircraftProperties.WingboxMaterial["density"]
    __shear = AircraftProperties.WingboxMaterial["shear strength"]

    integrationLimit = 50
    integrationFidelity = 20

    #spar flange connection stringer should be the name basically with the orientation such aconnection stringer would have in the top right
    #stringers should be alist of lists = [[absstarty, absendy, relposx, stringerShape], [......], ...]
    #flange thickness should be list = [[[topt, bott], endposY], [[topt, bott]], endposY], .....], so first index starts at wingroot
    def __init__(self, ribLocations, spars, stringersTop, stringersBottom, sparFlangeConnectionStringerShape, flangeThicknesses, crosssectionAmount, wingLoading=WingLoads([], [None, None, None], [None, None, None])):
        self.span = AircraftProperties.Planform["span"]
        self.semispan = self.span/2
        self.taperratio = AircraftProperties.Planform["taper ratio"]
        self.rootchord = AircraftProperties.Planform["root chord"]
        self.tipchord = self.rootchord * self.taperratio
        self.leadingedgeXPosTip = self.leadingEdgeXPos(self.semispan) # x position of the leading edge of the tip chord

        #ribs
        self.ribLocations = ribLocations
        self.ribLines = []
        for rib in self.ribLocations:
            if rib < self.semispan:
                lexpos = self.leadingEdgeXPos(rib)
                xVals = [lexpos, lexpos + self.getChordAtY(rib)]
                yVals = [rib, rib]
                self.ribLines.append([xVals, yVals])

        #spars
        self.spars = spars
        self.sparLines = []
        for spar in self.spars:
            spar[1] = min(spar[1], self.semispan)
            xVals = [spar[2] * self.getChordAtY(spar[0]) + self.leadingEdgeXPos(spar[0]), spar[2] * self.getChordAtY(spar[1]) + self.leadingEdgeXPos(spar[1])]
            yVals = [spar[0], spar[1]]
            self.sparLines.append([xVals, yVals])

        #stringers
        self.stringersTop = stringersTop
        self.stringersBottom = stringersBottom
        self.stringers = [stringersTop, stringersBottom]
        self.connectionStringers = self.__generateSparFlangeConnectionStringers(sparFlangeConnectionStringerShape)

        self.stringerLines = []
        for stringerList in self.stringers:
            stringerLineList = []
            for stringer in stringerList:
                stringer[1] = min(stringer[1], self.semispan)
                xVals = [self.leadingEdgeXPos(stringer[0]) + stringer[2] * self.getChordAtY(stringer[0]), self.leadingEdgeXPos(stringer[1]) + stringer[2] * self.getChordAtY(stringer[1])]
                yVals = [stringer[0], stringer[1]]
                stringerLineList.append([xVals, yVals, stringer[3]])
            self.stringerLines.append(stringerLineList)

        for index, stringerList in enumerate(self.connectionStringers):
            for stringer in stringerList:
                xVals = [stringer[2], stringer[3]]
                yVals = [stringer[0], stringer[1]]
                self.stringerLines[index].append([xVals, yVals, stringer[4]])

        self.flangeThicknesses = flangeThicknesses

        self.planformPolygon = Polygon([[0, 0], [self.rootchord, 0], [self.leadingedgeXPosTip + self.tipchord, self.semispan], [self.leadingedgeXPosTip, self.semispan]])

        self.crosssectionAmount = crosssectionAmount

        self.crossectionYLocations = np.linspace(0, self.semispan, self.crosssectionAmount, endpoint=True)
        self.crosssecctions = self.generateCrossections()

        #list = [yLocations, Inertia]
        self.ixxList, self.izzList, self.ixzList, self.jList = self.__generateInertialists()

        self.ixxFuncY = interp1d(self.ixxList[0], self.ixxList[1])
        self.izzFuncY = interp1d(self.izzList[0], self.izzList[1])
        self.ixzFuncY = interp1d(self.ixzList[0], self.ixzList[1])
        self.jFuncY = interp1d(self.jList[0], self.jList[1])

        #just for drawing the inertias
        self.inertiaFunctionsY = [self.ixxFuncY, self.izzFuncY, self.ixzFuncY, self.jFuncY]

        #mat properties
        # shear modulus": 26 * 10**9,
        self.E = AircraftProperties.WingboxMaterial["e modulus"]
        self.G = AircraftProperties.WingboxMaterial["shear modulus"]

        self.centroidDistFromC4AtYFunc = self.centroidDistFromC4AtY()
        self.internalVolume, self.materialVolume = self.__getInternalAndMaterialVolume()

        self.totalMass = self.materialVolume * self.__density

        self.wingLoading = wingLoading

    def __getInternalAndMaterialVolume(self):
        internal = 0
        material = 0
        ydist = self.crossectionYLocations[1]-self.crossectionYLocations[0]

        for crosssection in self.crosssecctions:
            internal += crosssection.internalArea
            material += crosssection.materialArea

        internal *= ydist
        material *= ydist

        return internal, material

    def getStructuralMassDistribution(self):
        structuralMass = []

        for crosssection in self.crosssecctions:
            structuralMass.append(crosssection.materialArea * self.__density)

        return interp1d(self.crossectionYLocations, structuralMass)

    def getFuelMassDistribution(self, evenlyDistributed=True):
        fuelMass = []
        totalVolumeSoFar = 0
        ydist = self.crossectionYLocations[1]-self.crossectionYLocations[0] # assumes equal crosssection spacing

        for crosssection in self.crosssecctions:
            totalVolumeSoFar += ydist * crosssection.internalArea

            if totalVolumeSoFar < self.__fuelVolumeReq or evenlyDistributed:
                fuelMass.append(crosssection.internalArea * self.__fuelDensity)
            else:
                fuelMass.append(0)

        if evenlyDistributed:
            factor = (self.__fuelVolumeReq/2)/totalVolumeSoFar

            fuelMass = [i * factor for i in fuelMass]

        return interp1d(self.crossectionYLocations, fuelMass)

    def centroidDistFromC4AtY(self):
        return interp1d(self.crossectionYLocations, [i.centroid[0]-i.chordLength*0.25 for i in self.crosssecctions])

    def __generateInertialists(self):
        yLocation = []
        ixx = []
        izz = []
        ixz = []
        j = []

        for crosssection in self.crosssecctions:
            yLocation.append(crosssection.yLocation)

            ixx.append(crosssection.ixx)
            izz.append(crosssection.izz)
            ixz.append(crosssection.izx)
            j.append(crosssection.j)

        return [yLocation, ixx], [yLocation, izz], [yLocation, ixz], [yLocation, j]

    def __generateSparFlangeConnectionStringers(self, stringerType):
        stringersTop = []
        stringersBottom = []

        leftSideStringer = [stringerType, stringerType.getMirrorStringerX()]
        rightSideStringer = [stringerType.getMirrorStringerZ(), stringerType.getMirrorStringerZ().getMirrorStringerX()]
        stringerPerSide = [leftSideStringer, rightSideStringer]

        directions = [-1, 1]

        for index, sparLine in enumerate(self.sparLines):
            sides = [0, 1]
            factor = 0.5

            if index == 0: # if left spar
                factor = 1
                sides = [1]
            elif index == 1: # if right spar
                factor = 1
                sides = [0]

            for side in sides:
                startx = sparLine[0][0] + self.spars[index][3](sparLine[1][0]) * directions[side] * factor + stringerType.baseLength/2 * directions[side]
                endx = sparLine[0][1] + self.spars[index][3](sparLine[1][1]) * directions[side] * factor + stringerType.baseLength/2 * directions[side]

                stringersTop.append([sparLine[1][0], sparLine[1][1], startx, endx, stringerPerSide[side][0]])
                stringersBottom.append([sparLine[1][0], sparLine[1][1], startx, endx, stringerPerSide[side][1]])
        return [stringersTop, stringersBottom]

    def generateCrossections(self):
        crosssections = []

        for yPos in self.crossectionYLocations:
            crosssections.append(self.getCrosssectionAtY(yPos))

        return crosssections

    def getGeneratedCrosssectionAtY(self, yPos):
        i = int((yPos * self.crosssectionAmount)/self.semispan)
        return self.crosssecctions[i]

    def getCrosssectionAtY(self, posY):
        flangeThicknesses = [self.flangeThicknesses[0](posY), self.flangeThicknesses[1](posY)]

        stringers = [] #top, bottom
        for stringerLineList in self.stringerLines:
            stringerPerSide = []
            for stringer in stringerLineList:
                xPos = self.__getLineIntersection(stringer[0], stringer[1], posY)
                if xPos is not None:
                    stringerPerSide.append([xPos - self.leadingEdgeXPos(posY), stringer[2]])
            stringers.append(stringerPerSide)

        sparThicknesses = []
        spars = []
        for index, spar in enumerate(self.sparLines):
            xPos = self.__getLineIntersection(spar[0], spar[1], posY)
            if xPos is not None:
                spars.append(xPos - self.leadingEdgeXPos(posY))
                sparThicknesses.append(self.spars[index][3](posY))

        return WingboxCrossection(chordLength=self.getChordAtY(posY), sparLocations=spars, sparThicknesses=sparThicknesses,
                                  flangeThicknesses=flangeThicknesses, stringersTop=stringers[0], stringersBottom=stringers[1], yLocation=posY)

    def __getLineIntersection(self, xList, yList, posY):
        if xList[0] == xList[1]: # this happens if at the center line of the wing, ugly way to do it but it works so w/e
            return xList[0]

        if yList[0] <= posY <= yList[1]:
            k = (yList[0] - yList[1]) / (xList[0] - xList[1])
            d = yList[0] - k * xList[0]

            return (posY-d)/k
        else:
            return None

    def getChordAtY(self, y):
        return self.rootchord * (1 + ((2 * (self.taperratio - 1)) / self.span) * y)

    def leadingEdgeXPos(self, posY):
        return self.rootchord / 2 - self.getChordAtY(posY) / 2

    def getVerticalDeflectionFunction(self, fidelity=integrationFidelity, limit=integrationLimit):
        MxAtYFunction = self.wingLoading.getInternalMoment(0)

        dvoverdy = []
        v = []
        yList = np.linspace(0, self.semispan, fidelity)

        for yVal in yList:
            f = lambda y: MxAtYFunction(y) / (self.E * self.ixxFuncY(y))
            i = quad(f, 0, yVal, limit=limit)
            dvoverdy.append(i[0] * -1)

        dvoverdyfunc = interp1d(yList, dvoverdy, kind='cubic')

        for yVal in np.linspace(0, self.semispan, fidelity):
            i = quad(dvoverdyfunc, 0, yVal, limit=limit)
            v.append(i[0])

        return interp1d(yList, v, kind='cubic')

    def getVerticalDeflectionAtY(self, yPos, fidelity=integrationFidelity, limit=integrationLimit):
        MxAtYFunction = self.wingLoading.getInternalMoment(0)

        dvoverdy = []
        yList = np.linspace(0, self.semispan, fidelity)

        for yVal in yList:
            f = lambda y: MxAtYFunction(y) / (self.E * self.ixxFuncY(y))
            i = quad(f, 0, yVal, limit=limit)
            dvoverdy.append(i[0] * -1)

        dvoverdyfunc = interp1d(yList, dvoverdy, kind='cubic')

        i = quad(dvoverdyfunc, 0, yPos, limit=limit)

        return i[0]

    def getTwistFunction(self, fidelity=integrationFidelity, limit=integrationLimit):
        MyAtYFunction = self.wingLoading.getInternalMoment(1)

        dthetaoverdy = []
        theta = []
        yList = np.linspace(0, self.semispan, fidelity)

        for yVal in yList:
            f = lambda y: MyAtYFunction(y) / (self.G * self.jFuncY(y))
            i = quad(f, 0, yVal, limit=limit)
            dthetaoverdy.append(i[0] * -1)

        dthethaoverdyfunc = interp1d(yList, dthetaoverdy, kind='cubic')

        for yVal in np.linspace(0, self.semispan, fidelity):
            i = quad(dthethaoverdyfunc, 0, yVal, limit=limit)
            theta.append(i[0])

        return interp1d(yList, theta, kind='cubic')

    def getTwist(self, MyAtYFunction, yPos, fidelity=integrationFidelity, limit=integrationLimit):
        MyAtYFunction = self.wingLoading.getInternalMoment(1)

        dthetaoverdy = []
        yList = np.linspace(0, self.semispan, fidelity)

        for yVal in yList:
            f = lambda y: MyAtYFunction(y) / (self.G * self.jFuncY(y))
            i = quad(f, 0, yVal, limit=limit)
            dthetaoverdy.append(i[0] * -1)

        dthethaoverdyfunc = interp1d(yList, dthetaoverdy, kind='cubic')

        i = quad(dthethaoverdyfunc, 0, yPos, limit=limit)

        return i[0]

    def getMaximumMinimumNormalStress(self):
        min = [0, None, [0, 0]] # magnitude, crosssection
        max = [0, None, [0, 0]] # magnitude, crosssection

        for crosssection in self.crosssecctions:
            for point in crosssection.outsidePolygon.coords:
                stress = crosssection.getBendingStressAtPoint(self.wingLoading.getInternalMoment(0)(crosssection.yLocation), self.wingLoading.getInternalMoment(2)(crosssection.yLocation), point[0], point[1])
                if stress > max[0]:
                    max = [stress, crosssection, point]
                elif stress < min[0]:
                    min = [stress, crosssection, point]

        return max, min

    def getMaximumShearStressMagnitude(self):
        max = -1 # magnitude
        yLoc = -1

        for crosssection in self.crosssecctions:
            stress = crosssection.getMaxShearStress(self.wingLoading.getShearForce(2)(crosssection.yLocation), self.wingLoading.getInternalMoment(1)(crosssection.yLocation))
            if stress > max:
                max = stress
                yLoc = crosssection.yLocation

        return max, yLoc

    def draw(self, drawTopStringers=True, drawBottomStringers=True):
        plt.clf()


        #Planform
        rotatedPlanform = self.planformPolygon.getRotatedPolygon(np.deg2rad(90))
        rotatedPlanform = rotatedPlanform.getTranslatedPolygon([-rotatedPlanform.coords[1][i] for i in range(2)])
        rotatedPlanform.addToPlot(plt, color="red")

        #draw ribs
        for ribline in self.ribLines:
            plt.plot(ribline[1], ribline[0], color="black")

        #draw spars
        for sparLine in self.sparLines:
            plt.plot(sparLine[1], sparLine[0], color="blue")

        #draw stringers
        if drawTopStringers:
            for stringer in self.stringerLines[0]:
                plt.plot(stringer[1], stringer[0], color="green")

        if drawBottomStringers:
            for stringer in self.stringerLines[1]:
                plt.plot(stringer[1], stringer[0], color="pink")

        # x, y ranges and same scale
        plt.ylim(-0.05*self.rootchord, 1.05 * self.rootchord)
        plt.xlim(-0.05 * self.semispan, 1.05 * self.semispan)
        plt.gca().set_aspect('equal', adjustable='box')

        plt.show()

    def drawCrosssection(self, posY, drawSidewallCenterlines=False, drawCentroid=False, drawBendingStress=False):
        self.getGeneratedCrosssectionAtY(posY).drawCrosssection(drawSidewallCenterlines, drawCentroid, drawBendingStress, Mx=self.wingLoading.getInternalMoment(0)(posY), Mz=self.wingLoading.getInternalMoment(2)(posY))

    def drawInertias(self):
        inertiaNames = [r"$I_{xx}$", r"$I_{zz}$", r"$I_{xz}$", "$J$"]
        fig, plots = plt.subplots(2, 2)
        fig.suptitle('Inertias')
        for row in range(2):
            for col in range(2):
                plots[row, col].plot(self.ixxList[0], [self.inertiaFunctionsY[row * 2 + col](i) for i in self.ixxList[0]])
                plots[row, col].set_title(inertiaNames[row * 2 + col])
                plots[row, col].set(ylabel=r'Inertia [$m^4$]', xlabel='semi-span [m]')
        fig.tight_layout(pad=1.5)

        plt.show()

    def drawDeflection(self, verticalDeflectionFunction, axisSameScale=False):
        plt.title("Vertical Deflection")
        plt.plot(self.ixxList[0], [verticalDeflectionFunction(i) for i in self.ixxList[0]])
        plt.xlabel('semi-span [m]')
        plt.ylabel('deflection [m]')

        if axisSameScale:
            plt.gca().set_aspect('equal', adjustable='box')

        plt.show()

    def drawTwist(self, twistDeflectionFunction, degrees=False):
        plt.title("Twist")
        if degrees:
            plt.plot(self.ixxList[0], [np.rad2deg(twistDeflectionFunction(i)) for i in self.ixxList[0]])
        else:
            plt.plot(self.ixxList[0], [twistDeflectionFunction(i) for i in self.ixxList[0]])
        plt.xlabel('semi-span [m]')
        if degrees:
            plt.ylabel('twist [deg]')
        else:
            plt.ylabel('twist [rad]')

        plt.show()

    def drawCentroidXatY(self):
        plt.title("Centroid X Position at Y")
        plt.plot(self.crossectionYLocations, [self.centroidDistFromC4AtYFunc(i) for i in self.crossectionYLocations])
        plt.xlabel('semi-span [m]')
        plt.ylabel(ylabel=r'x ccentroid Position [m]')

        plt.show()

    def drawStructuralMassDistribution(self):
        plt.title("Material Mass Distribution")
        structMassDistr = self.getStructuralMassDistribution()
        plt.plot(self.crossectionYLocations, [structMassDistr(i) for i in self.crossectionYLocations])
        plt.xlabel('semi-span [m]')
        plt.ylabel(ylabel=r'mass [$\frac{kg}{m}$]')

        plt.show()

    def drawFuelMassDistribution(self):
        plt.title("Fuel Mass Distribution")
        fuelMassDistr = self.getFuelMassDistribution()
        plt.plot(self.crossectionYLocations, [fuelMassDistr(i) for i in self.crossectionYLocations])

        plt.xlabel('semi-span [m]')
        plt.ylabel(ylabel=r'mass [$\frac{kg}{m}$]')

        plt.show()





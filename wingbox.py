import matplotlib.pyplot as plt
from aircraftProperties import AircraftProperties
from Polygon import Polygon
from wingboxCrosssection import WingboxCrossection
from scipy.interpolate import interp1d
from scipy.integrate import quad
from WingLoads import WingLoads
import numpy as np
import bucklingcoeff
import operator

class Wingbox:
    __fuelVolumeReq = AircraftProperties.Fuel["total fuel volume required"] / 2
    __fuelDensity = AircraftProperties.Fuel["fuel density"]
    __density = AircraftProperties.WingboxMaterial["density"]
    __shearStrength = AircraftProperties.WingboxMaterial["shear strength"]
    __yieldStrength = AircraftProperties.WingboxMaterial["yield strength"]
    __poissons = AircraftProperties.WingboxMaterial["poisson's ratio"]
    __E = AircraftProperties.WingboxMaterial["e modulus"]
    __G = AircraftProperties.WingboxMaterial["shear modulus"]


    integrationLimit = 50
    integrationFidelity = 20

    #spar flange connection stringer should be the name basically with the orientation such aconnection stringer would have in the top right
    #stringers should be alist of lists = [[absstarty, absendy, relposx, stringerShape], [......], ...]
    #flange thickness should be list = [[[topt, bott], endposY], [[topt, bott]], endposY], .....], so first index starts at wingroot
    def __init__(self, ribLocations, spars, stringersTop, stringersBottom, sparCapSide, sparCapCenter, flangeThicknesses, crosssectionAmount, ribThickness, wingLoading=WingLoads([], [None, None, None], [None, None, None])):
        self.span = AircraftProperties.Planform["span"]
        self.semispan = self.span/2
        self.taperratio = AircraftProperties.Planform["taper ratio"]
        self.rootchord = AircraftProperties.Planform["root chord"]
        self.tipchord = self.rootchord * self.taperratio
        self.leadingedgeXPosTip = self.leadingEdgeXPos(self.semispan) # x position of the leading edge of the tip chord

        self.__ribThickness = ribThickness

        #ribs
        self.ribLocations = []
        for i in ribLocations:
            self.ribLocations.append(min(self.semispan, i))

        minThick = 100000
        for i in self.ribLocations:
            for j in range(2):
                if flangeThicknesses[j](i) < minThick:
                    minThick = flangeThicknesses[j](i)

        self.minFlangeThickness = minThick
        self.rivetHeadSize = self.minFlangeThickness * 0.9
        self.__minRivetPitch = self.rivetHeadSize * 4

        print(self.__minRivetPitch)

        self.ribLines = []
        for rib in self.ribLocations:
            if rib <= self.semispan:
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
        self.sparCaps = self.__placeSparCaps(sparCapSide, sparCapCenter)
        self.sideSparCaps = sparCapSide

        self.stringerLines = []
        for stringerList in self.stringers:
            stringerLineList = []
            for stringer in stringerList:
                stringer[1] = min(stringer[1], self.semispan)
                xVals = [self.leadingEdgeXPos(stringer[0]) + stringer[2] * self.getChordAtY(stringer[0]), self.leadingEdgeXPos(stringer[1]) + stringer[2] * self.getChordAtY(stringer[1])]
                yVals = [stringer[0], stringer[1]]
                stringerLineList.append([xVals, yVals, stringer[3]])
            self.stringerLines.append(stringerLineList)

        for index, stringerList in enumerate(self.sparCaps):
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

        self.centroidDistFromC4AtYFunc = self.centroidDistFromC4AtY()
        self.internalVolume, self.materialVolume = self.__getInternalAndMaterialVolume()

        self.totalMass = self.materialVolume * self.__density

        self.wingLoading = wingLoading

    #region buckling checks

    def checkShearWebBuckling(self):
        shearList = self.getMaximumShearStressList()

        maxShearPerSection = [[0] * (len(self.ribLocations) - 1), [0] * (len(self.ribLocations) - 1)] # list of lists [front, aft]

        aPerSection = []
        bPerSection = [[], []]

        sidewallIndices = [[0, 1], [2, 3]]

        for index, yLocation in enumerate(self.crossectionYLocations):
            sectionIndex = -1
            for secIndex, loc in enumerate(self.ribLocations):
                if yLocation >= loc:
                    sectionIndex = min(secIndex, len(self.ribLocations) - 2)


            #max shear magnitude
            for i in range(2):
                #shear
                if shearList[i][index] > maxShearPerSection[i][sectionIndex]:
                    maxShearPerSection[i][sectionIndex] = shearList[i][index]

        for sectionIndex in range(len(self.ribLocations) - 1):
            #a per section
            aPerSection.append(((self.ribLocations[sectionIndex + 1] - self.ribLocations[sectionIndex])**2 + (self.leadingEdgeXPos(self.ribLocations[sectionIndex + 1]) - self.leadingEdgeXPos(self.ribLocations[sectionIndex]))**2)**0.5)

            for i in range(2):
                # b, spar length at each end
                outerPolygon = self.getGeneratedCrosssectionAtY(self.ribLocations[sectionIndex]).outsidePolygon
                dist = abs(outerPolygon.coords[sidewallIndices[i][0]][1] - outerPolygon.coords[sidewallIndices[i][1]][1])
                bPerSection[i].append(dist)

        for sectionIndex in range(len(self.ribLocations)-1):

            sparThicknesses = [self.spars[0][3]((self.ribLocations[sectionIndex + 1] + self.ribLocations[sectionIndex])/2), self.spars[0][3]((self.ribLocations[sectionIndex + 1] + self.ribLocations[sectionIndex])/2)]

            for i in range(2):
                ks = bucklingcoeff.shearfuncC(aPerSection[sectionIndex]/bPerSection[i][sectionIndex])
                tao_cr = min(((np.pi**2 * ks * self.__E) / (12 * (1 - self.__poissons**2))) * (sparThicknesses[i] / bPerSection[i][sectionIndex]) ** 2, self.__shearStrength)
                if maxShearPerSection[i][sectionIndex] > tao_cr:
                    sideText = ["front", "aft"]
                    print("Shear Buckling occurs in section", str(sectionIndex + 1), "from root. In spar:", sideText[i], ". Max Shear stress in that sheet:", str(maxShearPerSection[i][sectionIndex]), "[Pa].  Maximum allowable shear stress:", str(tao_cr), "[Pa]")
                    return True

        return False

    def checkFlangeBuckling(self):
        #for checking if flange at top or obttom is in compression
        if self.wingLoading.getInternalMoment(0)(1) < 0:
            flangeThickness = self.flangeThicknesses[0]
            cornerCoords = [0, 3]
            stringerSideIndex = 0
            #print("Compression Top")
        else:
            flangeThickness = self.flangeThicknesses[1]
            cornerCoords = [1, 2]
            stringerSideIndex = 1
            #print("Compression Bot")

        compressiveStressList = self.getMaximumCompressiveStressMagnitudeList()

        maxCompForcePerSection = [0] * (len(self.ribLocations) - 1)

        for index, yLocation in enumerate(self.crossectionYLocations):
            sectionIndex = -1
            for secIndex, loc in enumerate(self.ribLocations):
                if yLocation >= loc:
                    sectionIndex = min(secIndex, len(self.ribLocations) - 2)

            #get flange area
            wingboxEdgesPolygon = self.getGeneratedCrosssectionAtY(yLocation).outsidePolygon

            flangeArea = ((wingboxEdgesPolygon.coords[cornerCoords[0]][0] - wingboxEdgesPolygon.coords[cornerCoords[1]][0])**2 + (wingboxEdgesPolygon.coords[cornerCoords[0]][1] - wingboxEdgesPolygon.coords[cornerCoords[1]][1])**2)**0.5 * flangeThickness(yLocation)

            force = compressiveStressList[index] * flangeArea
            if force > maxCompForcePerSection[sectionIndex]:
                maxCompForcePerSection[sectionIndex] = force

        for sectionIndex in range(len(self.ribLocations) - 1):
            #a per section
            a = self.ribLocations[sectionIndex + 1] - self.ribLocations[sectionIndex]

            crossection = self.getGeneratedCrosssectionAtY(self.ribLocations[sectionIndex])

            stringerPolygons = crossection.stringerPolygons[stringerSideIndex]
            stringerTypes = crossection.stringers[stringerSideIndex]
            skipNext = False
            for index in range(1, len(stringerPolygons)):
                if skipNext:
                    skipNext = False
                    continue

                indeces = [index, index - 1]
                clampingCounter = 0

                # check if index is sparcap, if so skip next to true
                if stringerTypes[index][1].isSparCap:
                    skipNext = True

                # check if one side is clamping, add to clamping counter
                for i in indeces:
                    if stringerTypes[i][1].isClampedAttachment:
                        clampingCounter += 1

                if clampingCounter == 0: # all sides simply supported
                    kcFunc = bucklingcoeff.compressionfuncH
                elif clampingCounter == 1: # one side simply supported, one clamped
                    kcFunc = bucklingcoeff.compressionfuncC
                else: #voth sides clamped
                    kcFunc = bucklingcoeff.compressionfuncCC

                t = flangeThickness(self.ribLocations[sectionIndex + 1])

                # b
                b = self.__getDistanceBetweenPoints(stringerPolygons[indeces[1]].referencePoints[3], stringerPolygons[indeces[0]].referencePoints[3])

                #checl if buckling
                kc = kcFunc(a / b)
                F_cr = ((np.pi**2 * kc * self.__E) / (12 * (1 - self.__poissons**2))) * (t / b) ** 2

                if maxCompForcePerSection[sectionIndex] > F_cr:
                    print("Skin Buckling occurs in section", str(sectionIndex + 1), "from root. Between stringers (including spar caps, starting at 1):", str(index), "and", str(index + 1),
                          "from left. Distance", str(b), "[m]. Max comp Force in that sheet:", str(maxCompForcePerSection[sectionIndex]),
                          "[N].  Maximum allowable Force:", str(F_cr), "[N]")
                    return True
        return False

    def checkColumnBuckling(self):
        K = 1 # both ends are pinned

        #for checking if flange at top or obttom is in compression
        if self.wingLoading.getInternalMoment(0)(1) < 0:
            stringerSideIndex = 0
        else:
            stringerSideIndex = 1

        maxInternalMomentPerSection = [[0, 0] for i in range(len(self.ribLocations))] # [yLoc, Mx]

        minReqPitch = [self.semispan, -1] # [pitch, sectionId]

        for index, yLocation in enumerate(self.crossectionYLocations):
            sectionIndex = -1
            for secIndex, loc in enumerate(self.ribLocations):
                if yLocation >= loc:
                    sectionIndex = min(secIndex, len(self.ribLocations) - 2)

            internalMoment = self.wingLoading.getInternalMoment(0)(yLocation)
            if abs(internalMoment) > abs(maxInternalMomentPerSection[sectionIndex][1]):
                maxInternalMomentPerSection[sectionIndex] = [yLocation, internalMoment]

        for sectionIndex in range(len(self.ribLocations)-1):
            crossection = self.getGeneratedCrosssectionAtY(maxInternalMomentPerSection[sectionIndex][0])

            stringerPolygons = crossection.stringerPolygons[stringerSideIndex]

            for index, stringer in enumerate(stringerPolygons):
                for j in range(2):
                    stress = crossection.getBendingStressAtPoint(Mx=maxInternalMomentPerSection[sectionIndex][1], Mz=0, x=stringer.referencePoints[j][0], z=stringer.referencePoints[j][1])
                    reqPitch = ((K * np.pi**2 * self.__E * stringer.getIxx())/(abs(stress)))**0.5
                    if reqPitch < minReqPitch[0]:
                        minReqPitch = [reqPitch, sectionIndex]

        if minReqPitch[0] < self.__minRivetPitch:
            print("Minimumm required rivet pitch of", str(minReqPitch[0]*1000), "[mm] is lower than minimum allowable pitch of", str(self.__minRivetPitch), "[mm] in section", str(minReqPitch[1] + 1), ".")
            return True

        return False

    def checkInterRivetBuckling(self):
        if self.wingLoading.getInternalMoment(0)(1) < 0:
            flangeThickness = self.flangeThicknesses[0]
        else:
            flangeThickness = self.flangeThicknesses[1]

        compressiveStressList = self.getMaximumCompressiveStressMagnitudeList()

        maxCompStressPerSection = [[0, 0] for i in range(len(self.ribLocations))] # [stress, yLoc]

        for index, yLocation in enumerate(self.crossectionYLocations):
            sectionIndex = -1
            for secIndex, loc in enumerate(self.ribLocations):
                if yLocation >= loc:
                    sectionIndex = min(secIndex, len(self.ribLocations) - 2)

            stress = compressiveStressList[index]
            if stress > maxCompStressPerSection[sectionIndex][0]:
                maxCompStressPerSection[sectionIndex] = [stress, yLocation]

        for i in range(len(self.ribLocations) - 1):
            stress_cr = 0.9 * AircraftProperties.Rivets["c"] * self.__E * (flangeThickness(maxCompStressPerSection[i][1])/self.__minRivetPitch)**2
            if maxCompStressPerSection[i][0] >= stress_cr:
                print("Inter rivet buckling occurs in section", str(i), ". Allowable stress:", str(stress_cr), "[Pa]. Actual stress:", str(maxCompStressPerSection[i]), "[Pa]")
                return True

        return False

    #endregion

    def __getDistanceBetweenPoints(self, coord1, coord2):
        return ((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2) ** 0.5

    def __getInternalAndMaterialVolume(self):
        internal = 0
        material = 0
        ydist = self.crossectionYLocations[1]-self.crossectionYLocations[0]

        for crosssection in self.crosssecctions:
            internal += crosssection.internalArea
            material += crosssection.materialArea

        internal *= ydist
        material *= ydist

        #rib contributions
        for i in self.ribLocations:
            for j in self.getGeneratedCrosssectionAtY(i).insidePolygons:
                ribVolume = j.getArea() * self.__ribThickness

                internal -= ribVolume
                material += ribVolume

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

    def __placeSparCaps(self, sparCapsSide, sparCapCenter):
        stringersTop = []
        stringersBottom = []

        leftSideCap = [sparCapsSide, sparCapsSide.getMirrorStringerX()]
        rightSideCap = [sparCapsSide.getMirrorStringerZ(), sparCapsSide.getMirrorStringerZ().getMirrorStringerX()]
        sideCaps = [leftSideCap, rightSideCap]

        leftCenterCap = [sparCapCenter, sparCapCenter.getMirrorStringerX()]
        rightCenterCap = [sparCapCenter.getMirrorStringerZ(), sparCapCenter.getMirrorStringerZ().getMirrorStringerX()]
        centerCaps = [leftCenterCap, rightCenterCap]

        directions = [-1, 1]

        for index, sparLine in enumerate(self.sparLines):
            sides = [0, 1]
            factor = 0.5

            lengthOffset = sparCapCenter.baseLength / 2
            stringerPerSide = centerCaps
            if index == 0: # if left spar
                factor = 1
                sides = [1]
                lengthOffset = sparCapsSide.baseLength/2
                stringerPerSide = sideCaps
            elif index == 1: # if right spar
                factor = 1
                sides = [0]
                lengthOffset = sparCapsSide.baseLength/2
                stringerPerSide = sideCaps

            for side in sides:
                startx = sparLine[0][0] + self.spars[index][3](sparLine[1][0]) * directions[side] * factor + lengthOffset * directions[side]
                endx = sparLine[0][1] + self.spars[index][3](sparLine[1][1]) * directions[side] * factor + lengthOffset * directions[side]

                stringersTop.append([sparLine[1][0], sparLine[1][1], startx, endx, stringerPerSide[side][0]])
                stringersBottom.append([sparLine[1][0], sparLine[1][1], startx, endx, stringerPerSide[side][1]])
        return [stringersTop, stringersBottom]

    def generateCrossections(self):
        crosssections = []

        for yPos in self.crossectionYLocations:
            crosssections.append(self.getCrosssectionAtY(yPos))

        return crosssections

    def getGeneratedCrosssectionAtY(self, yPos):
        i = min(int((yPos * self.crosssectionAmount)/self.semispan), self.crosssectionAmount - 1)
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
        if yList[0] <= posY <= yList[1]:
            if xList[0] == xList[1]: # this happens if at the center line of the wing, ugly way to do it but it works so w/e
                return xList[0]
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
            f = lambda y: MxAtYFunction(y) / (self.__E * self.ixxFuncY(y))
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
            f = lambda y: MxAtYFunction(y) / (self.__E * self.ixxFuncY(y))
            i = quad(f, 0, yVal, limit=limit)
            dvoverdy.append(i[0] * -1)

        dvoverdyfunc = interp1d(yList, dvoverdy, kind='cubic')

        i = quad(dvoverdyfunc, 0, yPos, limit=limit)

        return i[0]

    def getTwistFunction(self, fidelity=integrationFidelity, limit=integrationLimit):
        MyAtYFunction = self.wingLoading.getInternalMoment(1)

        theta = []
        yList = np.linspace(0, self.semispan, fidelity)

        for yVal in yList:
            dthethaoverdyfunc = lambda y: MyAtYFunction(y) / (self.__G * self.jFuncY(y))
            i = quad(dthethaoverdyfunc, 0, yVal, limit=limit)
            theta.append(i[0] * -1)

        return interp1d(yList, theta, kind='cubic')

    def getTwist(self, yPos, limit=integrationLimit):
        MyAtYFunction = self.wingLoading.getInternalMoment(1)

        dthethaoverdyfunc = lambda y: MyAtYFunction(y) / (self.__G * self.jFuncY(y))
        i = quad(dthethaoverdyfunc, 0, yPos, limit=limit)

        return i[0]

    def getMaximumMinimumNormalStress(self):
        min = [0, None, [0, 0]] # magnitude, crosssection
        max = [0, None, [0, 0]] # magnitude, crosssection

        for crosssection in self.crosssecctions:
            for point in crosssection.outsidePolygon.coords:
                stress = crosssection.getBendingStressAtPoint(self.wingLoading.getInternalMoment(0)(crosssection.yLocation), 0, point[0], point[1])
                if stress > max[0]:
                    max = [stress, crosssection, point]
                elif stress < min[0]:
                    min = [stress, crosssection, point]

        return max, min

    def getMaximumTensileStressList(self):
        maxNormalStressList = []

        for crosssection in self.crosssecctions:
            max = 0
            for point in crosssection.outsidePolygon.coords:
                stress = crosssection.getBendingStressAtPoint(self.wingLoading.getInternalMoment(0)(crosssection.yLocation), 0, point[0], point[1])
                if stress > max:
                    max = stress

            maxNormalStressList.append(max)

        return maxNormalStressList

    def getMaximumCompressiveStressMagnitudeList(self):
        minNormalStressList = []

        for crosssection in self.crosssecctions:
            min = 0
            for point in crosssection.outsidePolygon.coords:
                stress = crosssection.getBendingStressAtPoint(self.wingLoading.getInternalMoment(0)(crosssection.yLocation), 0, point[0], point[1])
                if stress < min:
                    min = stress

            minNormalStressList.append(abs(min))

        return minNormalStressList

    def getMaximumShearStressList(self):
        front = []
        aft = []

        for crosssection in self.crosssecctions:
            stress = crosssection.getMaxShearStress(self.wingLoading.getShearForce(2)(crosssection.yLocation), self.wingLoading.getInternalMoment(1)(crosssection.yLocation))

            front.append(stress[0])
            aft.append(stress[1])

        return [front, aft]

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
        self.getGeneratedCrosssectionAtY(posY).drawCrosssection(drawSidewallCenterlines, drawCentroid, drawBendingStress, Mx=self.wingLoading.getInternalMoment(0)(posY), Mz=0)

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

    def drawMaximumTensileStress(self):
        plt.title("Maximum Tensile Stress")
        plt.plot(self.crossectionYLocations, self.getMaximumTensileStressList())
        plt.xlabel('semi-span [m]')
        plt.ylabel('stress [Pa]')

        plt.show()

    def drawMaximumNormalStress(self):
        plt.title("Maximum Normal Stress")
        plt.plot(self.crossectionYLocations, [max(self.getMaximumTensileStressList()[i], self.getMaximumCompressiveStressMagnitudeList()[i]) for i in self.crossectionYLocations])
        plt.xlabel('semi-span [m]')
        plt.ylabel('stress [Pa]')

        plt.show()

    def drawNormalStressSafetyMargin(self):
        plt.title("Normal Stress Safety Margin")
        plt.plot(self.crossectionYLocations, [self.__yieldStrength/max(self.getMaximumTensileStressList()[i], self.getMaximumCompressiveStressMagnitudeList()[i]) for i in self.crossectionYLocations])
        plt.xlabel('semi-span [m]')
        plt.ylabel('safety-margin [-]')

        plt.show()

    def drawEdgeCrackSafetyMargin(self):
        plt.title("Edge Crack Safety Margin")
        edgecrackStrength = AircraftProperties.WingboxMaterial["edge crack strength"]
        plt.plot(self.crossectionYLocations, [edgecrackStrength/i for i in self.getMaximumTensileStressList()])
        plt.xlabel('semi-span [m]')
        plt.ylabel('safety-margin [-]')

        plt.show()

    def drawCenterCrackSafetyMargin(self):
        plt.title("Center Crack Safety Margin")
        centercrackStrength = AircraftProperties.WingboxMaterial["center crack strength"]
        plt.plot(self.crossectionYLocations, [centercrackStrength/i for i in self.getMaximumTensileStressList()])
        plt.xlabel('semi-span [m]')
        plt.ylabel('safety-margin [-]')

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


    def drawMaxShearStress(self):
        yLabels = [r"$\tau_{front}$ [Pa]", r"$\tau_{aft}$ [Pa]"]
        data = self.getMaximumShearStressList()

        fig, plots = plt.subplots(1, 2)
        fig.suptitle("Max Shear Stress Magnitude")
        for col in range(2):
            plots[col].plot(self.crossectionYLocations, data[col])
            plots[col].set(ylabel=yLabels[col], xlabel='semi-span [m]')
        fig.tight_layout(pad=1.5)

        plt.show()




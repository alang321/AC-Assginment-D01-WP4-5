from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
from Polygon import Polygon
import airfoilDataParser
from aircraftProperties import AircraftProperties
import operator


class WingboxCrossection:
    __yield = AircraftProperties.WingboxMaterial["yield strength"]
    __ultimate = AircraftProperties.WingboxMaterial["ultimate strength"]

    def __init__(self, chordLength, sparLocations, sparThicknesses, flangeThicknesses, stringersTop, stringersBottom, yLocation=0):
        #spar locations
        self.chordLength = chordLength

        #y location of crossection, not used in calculation, just for reference if needed
        self.yLocation = yLocation

        #spars
        self.sparAbsoluteLocations = sorted(sparLocations) # m
        self.sparLocations = [i/self.chordLength for i in self.sparAbsoluteLocations]
        self.sparThicknesses = sorted(sparThicknesses)

        self.numSections = len(sparLocations) - 1 #number of Sections

        self.flangeThicknesses = flangeThicknesses

        self.airfoilData = airfoilDataParser.data

        #make a list with the sidewall thicknesses for each secction
        self.perSectionThicknesses = [] # 0 left, 1 bot, 2 right, 3 top
        for sectionIndex in range(self.numSections):
            sectionThicknesses = [0 for i in range(4)]
            #flanges
            sectionThicknesses[1] = self.flangeThicknesses[0]
            sectionThicknesses[3] = self.flangeThicknesses[1]
            #spars
            sectionThicknesses[0] = self.sparThicknesses[sectionIndex]
            sectionThicknesses[2] = self.sparThicknesses[sectionIndex + 1]

            self.perSectionThicknesses.append(sectionThicknesses)

        #number of stringers lists
        self.stringersTop = sorted(stringersTop, key=operator.itemgetter(0))
        self.stringersBottom = sorted(stringersBottom, key=operator.itemgetter(0))
        self.stringers = [self.stringersTop, self.stringersBottom]

        #outside corner Coordinates of the entire wingbox
        outerCornerCoordinates = self.__getOuterCornerCoordinates()
        self.outsidePolygon = Polygon(outerCornerCoordinates)

        #center line Polygons for each section
        centerLineCornerCoordinateList = self.__getInnerCornerCoordinates(outerCornerCoordinates, centerLine=True)
        self.centerLinePolygons = []
        for cornerCoordinates in centerLineCornerCoordinateList:
            self.centerLinePolygons.append(Polygon(cornerCoordinates))

        #inside Polygons for each section
        insideCornerCoordinateList = self.__getInnerCornerCoordinates(outerCornerCoordinates)
        self.insidePolygons = []
        for cornerCoordinates in insideCornerCoordinateList:
            self.insidePolygons.append(Polygon(cornerCoordinates))

        #Polygons for each stringer
        self.stringerPolygons = self.__placeStringers()

        #what values are added and substracted for centroid
        self.positiveCentroidContributions = []
        self.negativeCentroidContributions = []
        self.positiveCentroidContributions.append(self.outsidePolygon)
        self.positiveCentroidContributions.extend(self.stringerPolygons[0])
        self.positiveCentroidContributions.extend(self.stringerPolygons[1])
        self.negativeCentroidContributions.extend(self.insidePolygons)

        #calculate total enclosed area and crossectional area
        self.materialArea = self.__getMaterialArea()
        self.enclosedArea = self.__getEnclosedArea()
        self.internalArea = self.__getInternalArea()

        #get centroid
        self.centroid = self.__getCentroid()

        #what values are added and substracted for inertia
        self.positiveInertiaContributions = []
        self.negativeInertiaContributions = []
        self.onlySteinerTermContribution = []
        self.positiveInertiaContributions.append(self.outsidePolygon)
        self.positiveInertiaContributions.extend(self.stringerPolygons[0])
        self.positiveInertiaContributions.extend(self.stringerPolygons[1])
        self.negativeInertiaContributions.extend(self.insidePolygons)

        #calculate inertias
        self.ixx, self.izz = self.__getIxxIzz()
        self.izx = self.__getIzx()
        self.j = self.__getIyy() # calculates the iyy using the thin walled assumption, It therefore neglectes contributionof the stringers

    # uses the forward and aft spar as fraction of chord, returns the coordinates of the wing box corners as a fraction of the chord in a list of tuples,
    def __getOuterCornerCoordinates(self):
        coordinates = []
        outsideSpars = [self.sparLocations[0], self.sparLocations[-1]]

        for spar in outsideSpars:
            for intersection in self.__intersection(spar):
                coordinates.append([spar * self.chordLength, intersection * self.chordLength])
        coordinates[3], coordinates[2] = coordinates[2], coordinates[3] # switch some coords so they are in the correct order
        return coordinates

    def __getInnerCornerCoordinates(self, outerCornerCoordinates, centerLine=False):
        coordinateList = []

        topFlangAngle = self.__getAngleLineHorizontal(outerCornerCoordinates[0], outerCornerCoordinates[3])
        bottomFlangeAngle = self.__getAngleLineHorizontal(outerCornerCoordinates[1], outerCornerCoordinates[2])
        flangeAngles = [topFlangAngle, bottomFlangeAngle]
        flangeVertPosForwardSpar = [outerCornerCoordinates[0][1], outerCornerCoordinates[1][1]]

        for sectionId in range(self.numSections):
            sectionCoordinates = []
            for sparId in range(sectionId, sectionId+2):
                if centerLine is True:
                    horfactor = 0.5
                    vertfactor = -0.5
                else:
                    horfactor = 1.0
                    vertfactor = -1.0

                # check if spar is not outside spar, then substract 0.5 from the the horfactor, since the spar locations on the inside are at the center of the spar
                if not (sparId == self.numSections or sparId == 0):
                    horfactor -= 0.5

                # if at the right side of current section multiply factor by -1, so it is shifted to the left
                if sparId == (sectionId + 1):
                    horfactor *= -1

                for flangeId, intersectionVerticalValue in enumerate(self.__intersection(self.sparLocations[sparId])):
                    # flange factor negative if at bottom
                    if flangeId == 1:
                        vertfactor *= -1

                    # shift coordinates inside, taking into account angle of flanges
                    horshift = horfactor * self.sparThicknesses[sparId]
                    vertical = flangeVertPosForwardSpar[flangeId] + np.tan(flangeAngles[flangeId]) * (self.sparAbsoluteLocations[sparId] - self.sparAbsoluteLocations[0] + horshift) + vertfactor * self.flangeThicknesses[flangeId] # dont question it :/

                    sectionCoordinates.append([self.sparAbsoluteLocations[sparId] + horshift, vertical])

            sectionCoordinates[3], sectionCoordinates[2] = sectionCoordinates[2], sectionCoordinates[3] # switch some coords so they are in the correct order
            coordinateList.append(sectionCoordinates)
        return coordinateList

    # calculate the centroid and returns it ass coordinate
    def __getCentroid(self):
        centroid = [0, 0]

        for polygon in self.positiveCentroidContributions:
            polygonCentroid = polygon.getCentroid()
            for axis in range(2):
                centroid[axis] += polygonCentroid[axis] * polygon.getArea()

        for polygon in self.negativeCentroidContributions:
            polygonCentroid = polygon.getCentroid()
            for axis in range(2):
                centroid[axis] -= polygonCentroid[axis] * polygon.getArea()

        for axis in range(2):
            centroid[axis] /= self.materialArea

        return centroid

    # calculate the inertia
    def __getIxxIzz(self):
        totalInertia = [0, 0] #m^4

        for polygon in self.positiveInertiaContributions:
            steiner = polygon.getSteinerTermIxxIzz(self.centroid)
            polygonInertia = [polygon.getIxx(), polygon.getIzz()]
            for axis in range(2):
                totalInertia[axis] += steiner[axis]
                totalInertia[axis] += polygonInertia[axis]

        for polygon in self.negativeInertiaContributions:
            steiner = polygon.getSteinerTermIxxIzz(self.centroid)
            polygonInertia = [polygon.getIxx(), polygon.getIzz()]
            for axis in range(2):
                totalInertia[axis] -= steiner[axis]
                totalInertia[axis] -= polygonInertia[axis]

        for polygon in self.onlySteinerTermContribution:
            steiner = polygon.getSteinerTermIxxIzz(self.centroid)
            for axis in range(2):
                totalInertia[axis] += steiner[axis]

        return totalInertia

    def __getIzx(self):
        totalInertia = 0 #m^4

        for polygon in self.positiveInertiaContributions:
            totalInertia += polygon.getIxz()
            totalInertia += polygon.getSteinerTermIxz(self.centroid)

        for polygon in self.negativeInertiaContributions:
            totalInertia -= polygon.getIxz()
            totalInertia -= polygon.getSteinerTermIxz(self.centroid)

        for polygon in self.onlySteinerTermContribution:
            totalInertia += polygon.getSteinerTermIxz(self.centroid)

        return totalInertia

    #iyy, with multisection method from appendix d1.3 of the readder, the equations are put into a lin equation solver, q1, q2, ....qn, dtheta/dy
    def __getIyy(self):
        return self.__getMulticellTorsion(1)[-1]**-1

    def __getMulticellTorsion(self, Torque):
        totalInertia = 0 #m^4

        #left side
        a = []
        #right side
        b = [0 for i in range(self.numSections + 1)] # always results in 0 for
        b[0] = Torque # except for first row, where you apply torque of 1, see reader for explanation

        #add first row, different procedure
        firstrow = [0 for i in range(self.numSections + 1)]
        for index, polygon in enumerate(self.centerLinePolygons):
            firstrow[index] = polygon.getArea() * 2
        a.append(firstrow)

        #add another row for all additional sections
        for sectionIndex in range(self.numSections):
            row = [0 for i in range(self.numSections + 1)]

            firstterm = 1/(2 * self.centerLinePolygons[sectionIndex].getArea()) #vairable for the first term by which everything is multiplied
            lengthEachSide = []

            for pointIndex in range(4):
                nextPointIndex = (pointIndex + 1) % 4
                lengthEachSide.append(self.__getDistanceBetweenPoints(self.centerLinePolygons[sectionIndex].coords[pointIndex], self.centerLinePolygons[sectionIndex].coords[nextPointIndex])) #length of side

            #add all side terms with g and thickness for shear flow coefficient of current section
            sum = 0
            for pointIndex in range(4):
                length = lengthEachSide[pointIndex]
                thickness = self.perSectionThicknesses[sectionIndex][pointIndex]

                sum += length / thickness

            row[sectionIndex] = firstterm * sum # since the material is the same in each side, G is also applied here

            #add negative term to shear flow coefficient to q of sections in contact
            #check left
            if sectionIndex != 0: # if not leftmost section
                sideIndex = 0
                row[sectionIndex - 1] = -((lengthEachSide[sideIndex] / self.perSectionThicknesses[sectionIndex][sideIndex])*firstterm)
            #check right
            if sectionIndex != (self.numSections - 1): # if not rightmost section
                sideIndex = 2
                row[sectionIndex + 1] = -((lengthEachSide[sideIndex] / self.perSectionThicknesses[sectionIndex][sideIndex])*firstterm)

            row[-1] = -1
            a.append(row)
        x = np.linalg.solve(a, b)
        return x

    # calculate the centroid and returns it ass coordinate
    def __placeStringers(self):
        stringers = [[], []]

        attachmentLine = [[self.insidePolygons[0].coords[0], self.insidePolygons[0].coords[3]], [self.insidePolygons[0].coords[1], self.insidePolygons[0].coords[2]]]

        for side in range(2):
            xList = []
            for stringer in self.stringers[side]:
                xList.append(stringer[0])

            zList = self.__getZfromXLine(attachmentLine[side][0], attachmentLine[side][1], xList)

            for index, stringer in enumerate(self.stringers[side]):
                stringers[side].append(stringer[1].getStringerPolygonAtPointAndLine([xList[index], zList[index]], [attachmentLine[side][0], attachmentLine[side][1]]))

        return stringers

    # calculate the area enclosed by wingbox in terms of chord
    def __getEnclosedArea(self):
        sum = 0
        for i in self.centerLinePolygons:
            sum += i.getArea()
        return sum

    def __getInternalArea(self):
        sum = 0
        for i in self.insidePolygons:
            sum += i.getArea()
        return sum

    def __getMaterialArea(self):
        totalArea = 0

        for polygonPositive in self.positiveCentroidContributions:
            totalArea += polygonPositive.getArea()

        for polygonNegative in self.negativeCentroidContributions:
            totalArea -= polygonNegative.getArea()

        return totalArea

    #this method is not very general and assumes that the maximum shear stress occurs at the left spar on the neutral axis
    def getMaxShearStress(self, Vz, internalTorque, takeIntoAccountStringers=True):
        Q = 0

        outsideCut = self.outsidePolygon.getCutByY(self.centroid[1], getTop=True)
        #outsideCut.addToPlot(plt, drawCentroid=True, centroidText=False)
        Q += abs(outsideCut.getCentroid()[1] - self.centroid[1]) * outsideCut.getArea()

        for i in self.insidePolygons:
            cut = i.getCutByY(self.centroid[1], getTop=True)
            #cut.addToPlot(plt, drawCentroid=True, color="blue", centroidText=False)
            Q -= abs(cut.getCentroid()[1] - self.centroid[1]) * cut.getArea()

        if takeIntoAccountStringers:
            for i in self.stringerPolygons[0]:
                if i.getCentroid()[1] > self.centroid[1]:
                    Q += abs(i.getCentroid()[1] - self.centroid[1]) * i.getArea()
                    #i.addToPlot(plt, drawCentroid=True, centroidText=False)

        totalSparThickness = sum(self.sparThicknesses)

        taoShear = (Vz * Q) / (self.ixx * totalSparThickness)

        qs = self.__getMulticellTorsion(internalTorque)

        front = -1
        aft = -1

        aftSparIndex = self.numSections - 1

        # if torque is positive and shear is positive front spar
        if internalTorque >= 0 and Vz >= 0:
            front = abs(taoShear) + abs(qs[0] / self.sparThicknesses[0])
            aft = abs(taoShear) - abs(qs[aftSparIndex] / self.sparThicknesses[aftSparIndex])
        # if torque is negative and shear is positive at aft spar
        elif internalTorque <= 0 and Vz >= 0:
            front = abs(taoShear) - abs(qs[0] / self.sparThicknesses[0])
            aft = abs(taoShear) + abs(qs[aftSparIndex] / self.sparThicknesses[aftSparIndex])
        # if torque is negative and shear is negative at front spar
        elif internalTorque <= 0 and Vz <= 0:
            front = abs(taoShear) + abs(qs[0] / self.sparThicknesses[0])
            aft = abs(taoShear) - abs(qs[aftSparIndex] / self.sparThicknesses[aftSparIndex])
        # if torque is positive and shear is negative at aft spar
        elif internalTorque >= 0 and Vz <= 0:
            front = abs(taoShear) - abs(qs[0] / self.sparThicknesses[0])
            aft = abs(taoShear) + abs(qs[aftSparIndex] / self.sparThicknesses[aftSparIndex])

        #plt.gca().set_aspect('equal', adjustable='box')
        #plt.show()

        return [abs(front), abs(aft)]

    def getBendingStressAtPoint(self, Mx, Mz, x, z):
        return ((self.ixx * Mz - self.izx * Mx)/(self.ixx * self.izz - self.izx**2))*x + ((self.izz * Mx - self.izx * Mz)/(self.ixx * self.izz - self.izx**2))*z

    #return the intersection of both the upper surface and the lower surface at specified location
    def __intersection(self, location):
        f_upper = interp1d(self.airfoilData[0][0], self.airfoilData[0][1])
        f_lower = interp1d(self.airfoilData[1][0], self.airfoilData[1][1])

        return float(f_upper(location)), float(f_lower(location))

    #calculate the angle between a line defined by 2 coordinates and the x axis
    def __getAngleLineHorizontal(self, coord1, coord2):
        if coord2[0] == coord1[0]:
            return np.pi/2
        return np.arctan((coord2[1]-coord1[1])/(coord2[0]-coord1[0]))

    def __getDistanceBetweenPoints(self, coord1, coord2):
        return ((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2) ** 0.5

    def __getZfromXLine(self, coord1, coord2, xList):
        k = (coord1[1]-coord2[1])/(coord1[0]-coord2[0])
        d = coord1[1] - k * coord1[0]

        zList = []
        for x in xList:
            zList.append(k*x+d)

        return zList

    def addBendingStressHeatMapPlot(self, Mx, Mz, points=200, ultimateScale=False, yieldScale=False):
        xVals = []
        zVals = []
        intensity = []

        lines = [[self.outsidePolygon.coords[0], self.outsidePolygon.coords[3]],
                 [self.outsidePolygon.coords[1], self.outsidePolygon.coords[2]]]

        verticalLines = [[self.outsidePolygon.coords[0], self.outsidePolygon.coords[1]],
                         [self.outsidePolygon.coords[2], self.outsidePolygon.coords[3]]]

        if self.numSections > 1:
            for i in range(self.numSections-1):
                verticalLines.append([self.centerLinePolygons[i].coords[3], self.centerLinePolygons[i].coords[2]])

        for line in lines:
            #along top lines with certain x values
            xtemp = np.linspace(self.sparAbsoluteLocations[0], self.sparAbsoluteLocations[-1], points)
            xVals.extend(xtemp)
            zVals.extend(self.__getZfromXLine(line[0], line[1], xtemp))

        for line in verticalLines:
            #along side lines with certain z values
            for i in np.linspace(line[0][1], line[1][1], points):
                zVals.append(i)
                xVals.append(line[0][0])

        for index in range(len(xVals)):
            intensity.append(self.getBendingStressAtPoint(Mx, Mz, xVals[index], zVals[index]))

        if ultimateScale:
            norm = plt.Normalize(vmin=-self.__ultimate, vmax=self.__ultimate)
            plt.scatter(xVals, zVals, c=intensity, norm=norm, cmap="plasma")
        elif yieldScale:
            norm = plt.Normalize(vmin=-self.__yield, vmax=self.__yield)
            plt.scatter(xVals, zVals, c=intensity, norm=norm, cmap="plasma")
        else:
            plt.scatter(xVals, zVals, c=intensity, cmap="plasma")

        plt.colorbar()

    def drawBendingStressHeatmap(self, Mx, Mz, points=300, ultimateScale=False, yieldScale=False):
        self.addBendingStressHeatMapPlot(Mx, Mz, points, ultimateScale, yieldScale)
        plt.gca().set_aspect('equal', adjustable='box')

        plt.show()

    def drawCrosssection(self, drawSidewallCenterlines=False, drawCentroid=False, drawBendingStress=False, Mx=None, Mz=0):
        plt.figure()

        if drawBendingStress:
            if Mx is not None:
                self.addBendingStressHeatMapPlot(Mx, Mz, yieldScale=True)

        #plot top and bottom line
        plt.plot([i*self.chordLength for i in self.airfoilData[0][0]], [i*self.chordLength for i in self.airfoilData[0][1]], color="blue")
        plt.plot([i*self.chordLength for i in self.airfoilData[1][0]], [i*self.chordLength for i in self.airfoilData[1][1]], color="blue")
        plt.plot([self.airfoilData[0][0][-1]*self.chordLength, self.airfoilData[1][0][-1]*self.chordLength], [self.airfoilData[0][1][-1]*self.chordLength, self.airfoilData[1][1][-1]*self.chordLength], color="blue")

        #draw centroid
        if drawCentroid:
            plt.plot(self.centroid[0], self.centroid[1], marker='o', color="red", markersize=3)
            plt.annotate("Centroid", self.centroid, color="red")

        #plot wingbox
        for i in self.positiveInertiaContributions:
            i.addToPlot(plt, color="red")

        for i in self.negativeInertiaContributions:
            i.addToPlot(plt, color="blue")

        for i in self.onlySteinerTermContribution:
            i.addToPlot(plt, color="green")

        if drawSidewallCenterlines is True:
            for i in self.centerLinePolygons:
                i.addToPlot(plt, color="grey", linestyle="dashdot", linewidth=0.8)

        #x, y ranges and same scale
        plt.xlim(-0.005 * self.chordLength, 1.005 * self.chordLength)
        plt.ylim(-0.1 * self.chordLength, 0.1 * self.chordLength)
        plt.gca().set_aspect('equal', adjustable='box')

        plt.show()

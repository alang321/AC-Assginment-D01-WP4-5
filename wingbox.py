import matplotlib.pyplot as plt
from aircraftProperties import AircraftProperties
from Polygon import Polygon
from wingboxCrosssection import WingboxCrossection

class Wingbox:
    #spar flange connection stringer should be the name basically with the orientation such aconnection stringer would have in the top right
    #stringers should be alist of lists = [[absstarty, absendy, relposx, stringerShape], [......], ...]
    #flange thickness should be list = [[[topt, bott], endposY], [[topt, bott]], endposY], .....], so first index starts at wingroot
    # spar thickness list = [[[spar1t, spar2 t, ..], endposY], [....], ...]
    def __init__(self, ribLocations, sparLocations, sparThicknesses, stringersTop, stringersBottom, sparFlangeConnectionStringerShape, flangeThicknesses, crosssectionAmount):
        self.span = AircraftProperties.Planform["span"]
        self.semispan = self.span/2
        self.taperratio = AircraftProperties.Planform["taper ratio"]
        self.rootchord = AircraftProperties.Planform["root chord"]
        self.tipchord = self.rootchord * self.taperratio
        self.leadingedgeXPosTip = self.leadingEdgeXPos(self.semispan) # x position of the leading edge of the tip chord

        #ribs
        self.ribLocations = [i * self.semispan for i in ribLocations]
        self.ribLines = []
        for rib in self.ribLocations:
            lexpos = self.leadingEdgeXPos(rib)
            xVals = [lexpos, lexpos + self.chordAtY(rib)]
            yVals = [rib, rib]
            self.ribLines.append([xVals, yVals])

        #spars
        self.sparLocations = sparLocations
        self.sparThicknesses = sparThicknesses
        self.sparLines = []
        for spar in self.sparLocations:
            xVals = [spar * self.rootchord, self.leadingedgeXPosTip + spar*self.tipchord]
            yVals = [0, self.semispan]
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
                xVals = [self.leadingEdgeXPos(stringer[0]) + stringer[2] * self.chordAtY(stringer[0]), self.leadingEdgeXPos(stringer[1]) + stringer[2] * self.chordAtY(stringer[1])]
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

        self.crosssecctions = self.generateCrossection(self.crosssectionAmount)

        #list = [yLocations, Inertia]
        self.ixx, self.izz, self.ixz, self.iyy = self.__generateInertialists()

    def __generateInertialists(self):
        yLocation = []
        ixx = []
        izz = []
        ixz = []
        iyy = []

        for crosssection in self.crosssecctions:
            yLocation.append(crosssection.yLocation)

            ixx.append(crosssection.ixx)
            izz.append(crosssection.izz)
            ixz.append(crosssection.izx)
            iyy.append(crosssection.iyy)

        return [yLocation, ixx], [yLocation, izz], [yLocation, ixz], [yLocation, izz]

    def __generateSparFlangeConnectionStringers(self, stringerType):
        stringersTop = []
        stringersBottom = []

        leftSideStringer = [stringerType, stringerType.getMirrorStringerX()]
        rightSideStringer = [stringerType.getMirrorStringerZ(), stringerType.getMirrorStringerZ().getMirrorStringerX()]
        stringerPerSide = [leftSideStringer, rightSideStringer]

        directions = [-1, 1]

        for index, sparLocation in enumerate(self.sparLocations):
            sides = [0, 1]
            factor = 0.5

            if index == 0: # if right spar
                factor = 1
                sides = [1]
            elif index == len(self.sparLocations) - 1: # if left spar
                factor = 1
                sides = [0]

            for side in sides:
                startx = sparLocation * self.rootchord + self.sparThicknesses[0][0][index] * directions[side] * factor + stringerType.baseLength/2 * directions[side]
                endx = self.leadingedgeXPosTip + sparLocation * self.tipchord + self.sparThicknesses[-1][0][index] * directions[side] * factor + stringerType.baseLength/2 * directions[side]

                stringersTop.append([0, self.semispan, startx, endx, stringerPerSide[side][0]])
                stringersBottom.append([0, self.semispan, startx, endx, stringerPerSide[side][1]])
        return [stringersTop, stringersBottom]

    def generateCrossection(self, crosssectionAmount):
        crosssections = []

        for i in range(1, crosssectionAmount+1):
            yPos = self.semispan / crosssectionAmount * i
            crosssections.append(self.getCrosssectionAtY(yPos))

        return crosssections

    def getGeneratedCrosssectionAtY(self, yPos):
        i = int((yPos * self.crosssectionAmount)/self.semispan)
        return self.crosssecctions[i]

    def getCrosssectionAtY(self, posY):
        sparThicknesses = []
        for i in self.sparThicknesses:
            if i[1] >= posY:
                sparThicknesses = i[0]
                break

        flangeThicknesses = []
        for i in self.flangeThicknesses:
            if i[1] >= posY:
                flangeThicknesses = i[0]
                break

        stringers = [] #top, bottom
        for stringerLineList in self.stringerLines:
            stringerPerSide = []
            for stringer in stringerLineList:
                xPos = self.__getStringerIntersection(stringer, posY)
                if xPos is not None:
                    stringerPerSide.append([self.__getStringerIntersection(stringer, posY) - self.leadingEdgeXPos(posY), stringer[2]])
            stringers.append(stringerPerSide)

        return WingboxCrossection(chordLength=self.chordAtY(posY), sparLocations=self.sparLocations, sparThicknesses=sparThicknesses,
                           flangeThicknesses=flangeThicknesses, stringersTop=stringers[0], stringersBottom=stringers[1], yLocation=posY)



    def __getStringerIntersection(self, stringer, posY):
        xList = stringer[0]
        yList = stringer[1]
        if yList[0] <= posY <= yList[1]:
            k = (yList[0] - yList[1]) / (xList[0] - xList[1])
            d = yList[0] - k * xList[0]

            return (posY-d)/k
        else:
            return None


    def drawTopView(self):
        plt.clf()

        #Planform
        self.planformPolygon.addToPlot(plt, color="red")

        #draw ribs
        for ribline in self.ribLines:
            plt.plot(ribline[0], ribline[1], color="black")

        #draw spars
        for sparLine in self.sparLines:
            plt.plot(sparLine[0], sparLine[1], color="blue")

        #draw stringers
        for index, stringerLineList in enumerate(reversed(self.stringerLines)):
            if index == 1:
                color = "green"
            else:
                color = "pink"
            for stringer in stringerLineList:
                plt.plot(stringer[0], stringer[1], color=color)

        # x, y ranges and same scale
        plt.xlim(-0.1*self.rootchord, 1.1 * self.rootchord)
        plt.ylim(-0.1 * self.semispan, 1.1 * self.semispan)
        plt.gca().set_aspect('equal', adjustable='box')

        plt.show()

    def drawInertia(self, inertia, scatterPlot=False):
        plt.clf()

        if scatterPlot:
            plt.scatter(inertia[0], inertia[1], color="red")
        else:
            plt.plot(inertia[0], inertia[1], color="red")

        plt.show()

    def chordAtY(self, y):
        return self.rootchord * (1 + ((2 * (self.taperratio - 1)) / self.span) * y)

    def leadingEdgeXPos(self, posY):
        return self.rootchord/2 - self.chordAtY(posY)/2






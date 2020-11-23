import matplotlib.pyplot as plt
import numpy as np
import operator

class Polygon:
    #x horizontal, z vertical, y perpendicular
    def __init__(self, coordinateList, referencePoints=None):
        self.referencePoints = referencePoints
        if referencePoints is None:
            self.referencePoints = [[0, 0]]

        self.coords = coordinateList
        self.numCorners = len(self.coords)

        # for inertia calculation the coordinate list has to be in counter-clockwise order, this is checked by getting the signed Area
        signedArea = self.__calcSignedArea()
        if signedArea < 0:
            self.coords.reverse()

        self.__area = abs(signedArea)
        self.__centroid = None
        self.__ixx = None
        self.__izz = None
        self.__ixz = None

    #hide instance variable behind properties, so they are only calculated when actually called, and after the first call the value is saved
    #this is done so inertia, area and centroid don't have to be calculated on each rotation, translation or mirroring
    def getArea(self):
        return self.__area

    def getCentroid(self):
        if self.__centroid is None:
            self.__centroid = self.__calcCentroid()
        return self.__centroid

    def getIxx(self):
        if self.__ixx is None:
            self.__ixx = self.__calcInertia(1)
        return self.__ixx

    def getIzz(self):
        if self.__izz is None:
            self.__izz = self.__calcInertia(0)
        return self.__izz

    def getIxz(self):
        if self.__ixz is None:
            self.__ixz = self.__calcProductInertia()
        return self.__ixz

    def addToPlot(self, plt, color="red", linestyle="solid", linewidth=1.0, drawPointIndices=False, drawRefPoints=False):
        # plot point indices, should be in counter clockwise order for inertia to work
        if drawPointIndices:
            for index, point in enumerate(self.coords):
                plt.annotate(index, point)

        if drawRefPoints:
            for point in self.referencePoints:
                plt.plot(point[0], point[1], marker='o', color=color)

        coord = self.coords.copy()
        coord.append(coord[0])  # close polygon
        xlist, yList = zip(*coord)  # lists of x and y values
        plt.plot(xlist, yList, color=color, linestyle=linestyle, linewidth=linewidth)

    def draw(self, color="red", linestyle="solid", linewidth=1.0, drawPointIndices=False, drawRefPoints=False):

        plt.gca().set_aspect('equal', adjustable='box')
        self.addToPlot(plt, color=color, drawRefPoints=drawRefPoints, drawPointIndices=drawPointIndices, linestyle=linestyle, linewidth=linewidth)
        plt.show()

    def getMirroredPolygonX(self):
        transMatrix = [[1, 0],
                       [0, -1]]

        return self.getTransformedByMatrix(transMatrix)

    def getMirroredPolygonZ(self):
        transMatrix = [[-1, 0],
                       [0, 1]]

        return self.getTransformedByMatrix(transMatrix)

    def getRotatedPolygon(self, angle):
        transMatrix = [[np.cos(angle), -np.sin(angle)],
                       [np.sin(angle), np.cos(angle)]]

        return self.getTransformedByMatrix(transMatrix)

    def getTranslatedPolygon(self, vector):
        newCoords = []
        newRefs = []

        for coord in self.coords:
            newCoords.append(np.add(coord, vector))

        for point in self.referencePoints:
            newRefs.append(np.add(point, vector))

        return Polygon(newCoords, newRefs)

    def getTransformedByMatrix(self, transMatrix):
        newCoords = []
        newRefs = []

        for coord in self.coords:
            newCoords.append(np.dot(coord, transMatrix))

        for point in self.referencePoints:
            newRefs.append(np.dot(point, transMatrix))

        return Polygon(newCoords, newRefs)


    # dont question this absolute mess
    def getCutByY(self, y, getBottom=False, getTop=False):
        if getBottom == getTop:
            raise Exception("You can only get one side of the Polygon")

        newCoords = []

        if getBottom:
            comp = operator.lt
        else:
            comp = operator.ge

        i = -1

        for coord in range(self.numCorners):
            if comp(self.coords[coord][1], y):
                i = coord
                break

        if i == -1:
            return None

        intersectionsFound = 0
        initiali = i
        coordIndex = i % self.numCorners
        prevCoordIndex = (i - 1) % self.numCorners

        while True:

            if comp(self.coords[coordIndex][1], y):  # if on the desired side
                if intersectionsFound == 1:
                    coord1 = self.coords[coordIndex]
                    coord2 = self.coords[prevCoordIndex]
                    k = (coord1[0] - coord2[0]) / (coord1[1] - coord2[1])
                    d = coord1[0] - k * coord1[1]
                    newCoords.append([k * y + d, y])
                    intersectionsFound += 1
                    if coordIndex == initiali:
                        break
                newCoords.append(self.coords[coordIndex])
            elif intersectionsFound == 0:
                # this index and previous index
                coord1 = self.coords[coordIndex]
                coord2 = self.coords[prevCoordIndex]
                k = (coord1[0] - coord2[0]) / (coord1[1] - coord2[1])
                d = coord1[0] - k * coord1[1]
                newCoords.append([k * y + d, y])
                intersectionsFound += 1
            i += 1
            coordIndex = i % self.numCorners
            prevCoordIndex = (i - 1) % self.numCorners
            if coordIndex == initiali and intersectionsFound != 1:
                break
        return Polygon(newCoords)

    def copy(self):
        return Polygon(self.coords.copy(), self.referencePoints.copy())

    #steiner term for both x and z axis in a list
    def getSteinerTermIxxIzz(self, point):
        return [(self.getCentroid()[1] - point[1])**2 * self.getArea(), (self.getCentroid()[0] - point[0])**2 * self.getArea()]

    # steiner term for product inertia at point
    def getSteinerTermIxz(self, point):
        return (self.getCentroid()[0] - point[0]) * (self.getCentroid()[1] - point[1]) * self.getArea()

    def __calcCentroid(self):
        sumX = 0
        sumY = 0
        for index in range(self.numCorners):
            nextIndex = (index+1) % self.numCorners

            secondterm = (self.coords[index][0] * self.coords[nextIndex][1] - self.coords[nextIndex][0] * self.coords[index][1])
            sumX += (self.coords[index][0] + self.coords[nextIndex][0]) * secondterm
            sumY += (self.coords[index][1] + self.coords[nextIndex][1]) * secondterm

        return [sumX / (6*self.getArea()), sumY / (6*self.getArea())]

    def __calcSignedArea(self):
        sum = 0
        for index in range(self.numCorners):
            nextIndex = (index+1) % self.numCorners
            sum += self.coords[index][0] * self.coords[nextIndex][1] - self.coords[index][1] * self.coords[nextIndex][0]

        return sum/2

    # inertia around centroid http://richardson.eng.ua.edu/Former_Courses/CE_331_fa09/Projects/A_and_I_of_Polygon.pdf
    # 0 for z axis, 1 for x axis
    def __calcInertia(self, axisInverse):
        inertia = 0
        for index in range(self.numCorners):
            nextIndex = (index+1) % self.numCorners

            secondterm = (self.coords[index][0] * self.coords[nextIndex][1] - self.coords[nextIndex][0] *self.coords[index][1])

            inertia += (self.coords[index][axisInverse]**2 + self.coords[index][axisInverse]*self.coords[nextIndex][axisInverse]+self.coords[nextIndex][axisInverse]**2)*secondterm

        inertia /= 12

        #this method calculates around the reference axis, shift to centroidal parallel axis
        inertia -= self.getCentroid()[axisInverse]**2 * self.getArea()

        return inertia

    def __calcProductInertia(self):
        inertia = 0
        for i in range(self.numCorners):
            nexti = (i + 1) % self.numCorners

            secondterm = (self.coords[i][0] * self.coords[nexti][1] - self.coords[nexti][0] * self.coords[i][1])

            inertia += (self.coords[i][0] * self.coords[nexti][1] + 2 * self.coords[i][0] * self.coords[i][1] + 2 * self.coords[nexti][0] * self.coords[nexti][1]+self.coords[nexti][0]*self.coords[i][1])*secondterm

        inertia /= 24

        #this method calculates around the reference axis, shift to centroidal parallel axis
        inertia -= self.getCentroid()[0] * self.getCentroid()[1] * self.getArea()

        return inertia

# a stringer class with some extra properties to the polygon class, this class is created with a stringer shape and an attachment point, this is then sort of a blue print for further stringers called with the function
class StringerType:
    def __init__(self, coordinateList=None, baseAlignmentPoints=None, attachmentPoint=None, polygon=None): #attachmentLineIndices are the indices of the coordinates that form the "base line", so the line in contact with the sheet its attached to
        if polygon is not None:
            self.baseLength = self.__getDistanceBetweenPoints(polygon.referencePoints[0], polygon.referencePoints[1])
            self.stringerShape = polygon
        else:
            if len(baseAlignmentPoints) != 2 or baseAlignmentPoints[0] == baseAlignmentPoints[1]:
                raise Exception("Invalid Alignment Points")

            self.baseLength = self.__getDistanceBetweenPoints(baseAlignmentPoints[0], baseAlignmentPoints[1])

            if attachmentPoint is None: # if no attachment point is set, set the attachment point between the two alignment points
                attachmentPoint = [0, 0]
                for i in range(2):
                    attachmentPoint[i] = (baseAlignmentPoints[0][i] + baseAlignmentPoints[1][i]) / 2

            baseAlignmentPoints.append(attachmentPoint)
            self.stringerShape = Polygon(coordinateList, baseAlignmentPoints)

    #returns a polygon of the stringer with the attachemnt point at point and the alignemntpoints aligned rotationally with the line, lin eis supposed to be a list of 2 long lists (2 points)
    def getStringerPolygonAtPointAndLine(self, point, line, mirrorHorizontal=False, mirrorVertical=False):
        newStringerPolygon = self.stringerShape.copy()

        #mirror
        if mirrorHorizontal:
            newStringerPolygon = newStringerPolygon.getMirroredPolygonX()
        if mirrorVertical:
            newStringerPolygon = newStringerPolygon.getMirroredPolygonZ()

        #rotate to reference line
        newStringerPolygon = newStringerPolygon.getRotatedPolygon(self.__angleBetweenAlignmentAndLine(line))

        #translate to placement point
        newStringerPolygon = newStringerPolygon.getTranslatedPolygon([point[0] - newStringerPolygon.referencePoints[2][0], point[1] - newStringerPolygon.referencePoints[2][1]])

        return newStringerPolygon

    def drawUnplacedStringer(self, color="red"):
        self.stringerShape.draw(drawRefPoints=True, color=color)

    def getMirrorStringerX(self):
        return StringerType(polygon=self.stringerShape.getMirroredPolygonX())

    def getMirrorStringerZ(self):
        return StringerType(polygon=self.stringerShape.getMirroredPolygonZ())

    def __getDistanceBetweenPoints(self, coord1, coord2):
        return ((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2) ** 0.5

    def __angleBetweenAlignmentAndLine(self, line):
        refpoints = self.stringerShape.referencePoints
        angleLine = np.arctan((line[1][1] - line[0][1]) / (line[1][0] - line[0][0]))
        angleAlignment = np.arctan((refpoints[1][1] - refpoints[0][1]) / (refpoints[1][0] - refpoints[0][0]))

        return angleAlignment - angleLine

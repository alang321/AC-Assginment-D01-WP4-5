from aircraftProperties import AircraftProperties
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np


class WingboxCrossection:

    # sidewall thickness = [left, bottom, right, top], numStringers = [top, bottom]
    def __init__(self, chordLength, forwardSpar, aftSpar, numStringers, sidewallThickness, stringerArea):
        self.chordLength = chordLength
        self.numStringers = numStringers
        self.sidewallThickness = sidewallThickness

        self.stringerArea = stringerArea

        self.forwardSpar = forwardSpar
        self.aftSpar = aftSpar
        self.forwardSparAbsolute = forwardSpar * chordLength #m
        self.aftSparAbsolute = aftSpar * chordLength #m
        self.spars = [self.forwardSpar, self.aftSpar]
        self.sparsAbsolute = [self.forwardSparAbsolute, self.aftSparAbsolute] # m

        self.airfoilData = self.parseAirfoilData()

        self.cornerCoordinates = self.getCornerCoordinates()

        self.stringerLocations = self.getStringerLocations()

        #section Properties
        self.enclosedArea = self.getEnclosedArea()
        self.centroid = self.getCentroid()

        self.centroid = [0,0]

        self.ixx = self.getInertia() #horizontal
        self.izz = self.getInertia() #vertical
        self.iyy = self.getInertia() #polar

    # uses the forward and aft spar as fraction of chord, returns the coordinates of the wing box corners as a fraction of the chord in a list of tuples,
    def getCornerCoordinates(self):
        coordinates = []
        for spar in self.spars:
            for intersection in self.intersection(spar):
                coordinates.append([spar * self.chordLength, intersection * self.chordLength])
        coordinates[3], coordinates[2] = coordinates[2], coordinates[3] # switch some coords so they are in the correct order
        return coordinates

    # calculate the centroid and returns it ass coordinate
    def getCentroid(self):
        centroid = []
        # centroid of i stringers (area times distance)
        # centroid of trapezoid wingbox
        coords = self.getCornerCoordinates()
        # centroidx = coords[][] - coords[][]
        # centroidy = coords[][] - coords[][]

        # average of all centroids = final centroid
        return centroid

    # calculate the inertia
    def getInertia(self):
        inertia = 0 #m^4

        sidePanels = [] # centroid coords, angle in rad, length, thickness
        for index in range(4):
            nextIndex = (index+1)%4

            coordCentroid = [0, 0]
            for i in range(2):
                coordCentroid[i] = (self.cornerCoordinates[index][i] + self.cornerCoordinates[nextIndex][i]) / 2

            angle = WingboxCrossection.getAngleLineHorizontal(self.cornerCoordinates[index], self.cornerCoordinates[nextIndex])
            length = WingboxCrossection.getDistanceBetweenPoints(self.cornerCoordinates[index], self.cornerCoordinates[nextIndex])

            sidePanels.append([coordCentroid, angle, length, self.sidewallThickness[index]])
        print(sidePanels)

        return inertia

    # calculate the centroid and returns it ass coordinate
    def getStringerLocations(self):
        coordinateList = []
        coords = self.getCornerCoordinates()


        for numStringersPerSide in self.numStringers:
            x = coords[3][0]-coords[0][0]
            dx = x/(numStringersPerSide-1)

            for stringer in numStringersPerSide:
                dx * stringer

        
        return coordinateList

    # calculate the area enclosed by wingbox in terms of chord
    def getEnclosedArea(self):
        coords = self.getCornerCoordinates()
        fwdspar = coords[0][1] - coords[1][1]
        bckspar = coords[3][1] - coords[2][1]
        lengwingbox = coords[2][0] - coords[0][0]
        A_m = 0.5*lengwingbox*(fwdspar + bckspar)
        return A_m

    def drawWingbox(self, drawStringerLocations=False):
        #plot top and bottom line
        plt.plot([i*self.chordLength for i in self.airfoilData[0][0]], [i*self.chordLength for i in self.airfoilData[0][1]], color="blue")
        plt.plot([i*self.chordLength for i in self.airfoilData[1][0]], [i*self.chordLength for i in self.airfoilData[1][1]], color="blue")

        #plot wingbox
        coord = self.getCornerCoordinates()
        coord.append(coord[0])  # repeat the first point to create a 'closed loop'
        xs, ys = zip(*coord)  # create lists of x and y values
        plt.plot(xs, ys, color="red")

        #draw centroid
        plt.plot(self.centroid[0], self.centroid[1], marker='o', color='r', ls='')

        #x, y ranges and same scale
        plt.xlim(-0.1*self.chordLength, 1.1*self.chordLength)
        plt.ylim(-0.5*self.chordLength, 0.5*self.chordLength)
        plt.gca().set_aspect('equal', adjustable='box')

        plt.show()

    #return the intersection of both the upper surface and the lower surface at specified location
    def intersection(self, location):
        f_upper = interp1d(self.airfoilData[0][0], self.airfoilData[0][1])
        f_lower = interp1d(self.airfoilData[1][0], self.airfoilData[1][1])

        return float(f_upper(location)), float(f_lower(location))

    @staticmethod #calculate the angle between a line defined by 2 coordinates and the x axis
    def getAngleLineHorizontal(coord1, coord2):
        if coord2[0] == coord1[0]:
            return np.pi/2
        return np.arctan((coord2[1]-coord1[1])/(coord2[0]-coord1[0]))

    @staticmethod
    def getDistanceBetweenPoints(coord1, coord2):
        return ((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2) ** 0.5

    @staticmethod
    def parseAirfoilData():
        temp = []
        data = []

        file = "data/airfoil" + AircraftProperties.Airfoil["airfoil identifier"] + ".dat"
        with open(file) as f:
            lines = f.read().splitlines()

        emptyLineIndex = []
        #get index of empty lines
        for index, line in enumerate(lines):
            if line == "":
                emptyLineIndex.append(index)

        temp.append(lines[emptyLineIndex[0]+1:emptyLineIndex[1]])
        temp.append(lines[emptyLineIndex[1]+1:])

        for index, i in enumerate(temp):
            data.append([[], []])
            for line in i:
                splitline = line.lstrip().split()
                data[index][0].append(float(splitline[0]))
                data[index][1].append(float(splitline[1]))

        return data
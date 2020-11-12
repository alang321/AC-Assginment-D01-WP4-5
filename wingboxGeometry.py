import matplotlib

from aircraftProperties import AircraftProperties
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


class WingboxGeometry:

    def __init__(self, forwardSpar, aftSpar):
        self.forwardSpar = forwardSpar
        self.aftSpar = aftSpar
        self.spars = [self.forwardSpar, self.aftSpar]

        self.airfoilData = self.parseAirfoilData()

    # uses the forward and aft spar as fraction of chord, returns the coordinates of the wing box edges as a fraction of the chord in a list of tuples,
    def edgeCoordinates(self):
        coordinates = []
        for spar in self.spars:
            for intersection in self.intersection(spar):
                coordinates.append([spar, intersection])
        return coordinates

    def drawWingbox(self):
        #plot top and bottom line
        plt.plot(self.airfoilData[0][0], self.airfoilData[0][1], color="blue")
        plt.plot(self.airfoilData[1][0], self.airfoilData[1][1], color="blue")

        #plot wingbox
        coord = self.edgeCoordinates()
        coord[3], coord[2] = coord[2], coord[3] # switch some coords so they are in the correct order
        coord.append(coord[0])  # repeat the first point to create a 'closed loop'
        xs, ys = zip(*coord)  # create lists of x and y values
        plt.plot(xs, ys, color="red")

        #x, y ranges and same scale
        plt.xlim(-0.1, 1.1)
        plt.ylim(-0.5, 0.5)
        plt.gca().set_aspect('equal', adjustable='box')

        plt.show()

    #return the intersection of both the upper surface and the lower surface at specified location
    def intersection(self, location):
        f_upper = interp1d(self.airfoilData[0][0], self.airfoilData[0][1])
        f_lower = interp1d(self.airfoilData[1][0], self.airfoilData[1][1])

        return float(f_upper(location)), float(f_lower(location))

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




from aircraftProperties import AircraftProperties
from scipy.interpolate import interp1d

class WingboxGeometry:
    # takes the forward and aft spar as fraction of chord, returns the coordinates of the wing box edges as a fraction of the chord in a list of tuples,

    def __init__(self, forwardSpar, aftSpar):
        self.forwardSpar = forwardSpar
        self.aftSpar = aftSpar
        self.spars = [self.forwardSpar, self.aftSpar]

    def edgeCoordinates(self):
        coordinates = []
        for spar in self.spars:
            for intersection in WingboxGeometry.intersection(spar):
                coordinates.append((spar, intersection))
        return coordinates

    #return the intersection of both the upper surface and the lower surface at specified location
    @staticmethod
    def intersection(location):
        data = WingboxGeometry.parseAirfoilData()
        f_upper = interp1d(data[0][0], data[0][1])
        f_lower = interp1d(data[1][0], data[1][1])

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
            data.append([[],[]])
            for line in i:
                splitline = line.lstrip().split()
                data[index][0].append(float(splitline[0]))
                data[index][1].append(float(splitline[1]))

        return data




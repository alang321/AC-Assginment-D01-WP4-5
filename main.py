from aircraftProperties import AircraftProperties
from wingboxCrosssection import WingboxCrossection
import AerodynamicDataThingy

print(AircraftProperties.Fuselage["Fuselage length"])

wingbox = WingboxCrossection(forwardSpar=0.15, aftSpar=0.6, chordLength=10, numStringerTop=8, numStringersBottom=3, stringerInertia=0.000000002238636567, stringerArea=0.00005775, sidewallThickness=[0.1, 0.1, 0.1, 0.1])
print(wingbox.edgeCoordinates())
print(wingbox.calculateEnclArea())
wingbox.drawWingbox()

print(AerodynamicDataThingy.Dragacc(10))
AerodynamicDataThingy.drawgraphs()
from aircraftProperties import AircraftProperties
from wingboxCrosssection import WingboxCrossection
import AerodynamicDataThingy
from Polygon import Polygon
from Polygon import StringerType

print(AircraftProperties.Fuselage["Fuselage length"])

#stringer = StringerType([[0, 0], [0, -1/1000], [19/1000, -1/1000], [19/1000, -20/1000], [20/1000, -20/1000], [20/1000, 0]], [[0, 0], [20/1000, 0]])
stringer = StringerType([[0, 0], [0, -3/200], [17/200, -3/200], [17/200, -20/200], [20/200, -20/200], [20/200, 0]], [[0, 0], [20/200, 0]])
stringer.drawUnplacedStringer()

#2 section
#wingbox = WingboxCrossection(chordLength=10, sparLocations=[0.15, 0.3, 0.6], sparThicknesses=[0.1, 0.2, 0.1], flangeThicknesses=[0.05, 0.25], amountStringerTop=[9, 10], amountStringerBottom=[2, 3], stringerType=stringer)
#3 section
#wingbox = WingboxCrossection(chordLength=10, sparLocations=[0.15, 0.3, 0.5, 0.7], sparThicknesses=[0.1, 0.2, 0.15, 0.1], flangeThicknesses=[0.05, 0.25], amountStringerTop=[5, 5, 6], amountStringerBottom=[2, 3, 2], stringerType=stringer)
#single cell
wingbox = WingboxCrossection(chordLength=10, sparLocations=[0.15, 0.6], sparThicknesses=[0.01, 0.01], flangeThicknesses=[0.01, 0.01], amountStringerTop=[9], amountStringerBottom=[3], stringerType=stringer)


print("Enclosed Area:", wingbox.enclosedArea)
print("Crosssectional Area:", wingbox.totalCrossectionalArea)
print("Ixx:", wingbox.ixx)
print("Izz:", wingbox.izz)
print("Iyy:", wingbox.iyy)
print("Ixz:", wingbox.izx)
wingbox.drawWingbox(drawSidewallCenterlines=True, drawCentroid=True)




#print(AerodynamicDataThingy.Dragacc(10))
#AerodynamicDataThingy.drawgraphs()
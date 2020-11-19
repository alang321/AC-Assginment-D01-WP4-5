from aircraftProperties import AircraftProperties
from wingboxCrosssection import WingboxCrossection
import AerodynamicDataThingy
from Polygon import StringerType
from wingbox import Wingbox

print(AircraftProperties.Fuselage["Fuselage length"])

#stringer = StringerType([[0, 0], [0, -1/1000], [19/1000, -1/1000], [19/1000, -20/1000], [20/1000, -20/1000], [20/1000, 0]], [[0, 0], [20/1000, 0]])
stringer = StringerType([[0, 0], [0, -3/200], [17/200, -3/200], [17/200, -20/200], [20/200, -20/200], [20/200, 0]], [[0, 0], [20/200, 0]])
stringer.drawUnplacedStringer()

stringerTop = stringer
stringerBottom = stringer.getMirrorStringerX()


wingbox = Wingbox(ribLocations=[0, 0.1, 0.2, 0.3, 0.4, 0.5], sparLocations=[0.15, 0.6], sparThicknesses=[[[0.01, 0.01], 30]], stringersTop=[[0, 20, 0.3, stringer]], stringersBottom=[[0, 10, 0.4, stringerBottom]], sparFlangeConnectionStringerShape=stringer, flangeThicknesses=[[[0.01, 0.01], 30]], crosssectionAmount=400)
wingbox.drawTopView()

while True:
    pos = int(input("Crossection at location in m:"))
    print(wingbox.getGeneratedCrosssectionAtY(pos).yLocation)
    wingbox.getGeneratedCrosssectionAtY(pos).drawWingbox()
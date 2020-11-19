from aircraftProperties import AircraftProperties
from wingboxCrosssection import WingboxCrossection
from Polygon import StringerType
from wingbox import Wingbox

print(AircraftProperties.Fuselage["Fuselage length"])

#stringer = StringerType([[0, 0], [0, -1/1000], [19/1000, -1/1000], [19/1000, -20/1000], [20/1000, -20/1000], [20/1000, 0]], [[0, 0], [20/1000, 0]])
stringer = StringerType([[0, 0], [0, -3/200], [17/200, -3/200], [17/200, -20/200], [20/200, -20/200], [20/200, 0]], [[0, 0], [20/200, 0]])
stringer.drawUnplacedStringer()

stringerTop = stringer
stringerBottom = stringer.getMirrorStringerX()

# rib location : fraction location of semispan, spar location: fractional location of chord, sparThickness: [[[Thicknesses of spar], endpoint], [[Thicknesses of spar2], endpoint2]], stringersTop = [[starty, endy, fraction of chord, stringer]] flangeThickness = same as spar thickness
#wingbox = Wingbox(ribLocations=[0, 0.1, 0.2, 0.3, 0.4, 0.5], sparLocations=[0.15, 0.6], sparThicknesses=[[[0.01, 0.01], 30]], stringersTop=[[0, 20, 0.3, stringer]], stringersBottom=[[0, 10, 0.4, stringerBottom]], sparFlangeConnectionStringerShape=stringer, flangeThicknesses=[[[0.01, 0.01], 30]], crosssectionAmount=400)
wingbox = Wingbox(ribLocations=[0, 0.1, 0.2, 0.3, 0.4, 0.5], sparLocations=[0.15, 0.6], sparThicknesses=[[[0.01, 0.01], 30]], stringersTop=[], stringersBottom=[], sparFlangeConnectionStringerShape=stringer, flangeThicknesses=[[[0.01, 0.01], 30]], crosssectionAmount=400)
wingbox.drawTopView()

wingbox.drawInertias()

while True:
    pos = float(input("Crossection at location in m:"))
    print(wingbox.getCrosssectionAtY(pos).centroid[0])

    wingbox.getGeneratedCrosssectionAtY(pos).drawWingbox(drawCentroid=True, drawSidewallCenterlines=True)
from aircraftProperties import AircraftProperties
from wingboxGeometry import WingboxGeometry


print(AircraftProperties.Fuselage["Fuselage length"])

a = WingboxGeometry(0.15, 0.6)

print(a.edgeCoordinates())
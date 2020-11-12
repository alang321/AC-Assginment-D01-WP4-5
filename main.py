from aircraftProperties import AircraftProperties
from wingboxGeometry import WingboxGeometry


print(AircraftProperties.Fuselage["Fuselage length"])

wingbox = WingboxGeometry(forwardSpar=0.15, aftSpar=0.6)
test = WingboxGeometry(forwardSpar=0.087654, aftSpar=0.8)

print(wingbox.edgeCoordinates())

wingbox.drawWingbox()
test.drawWingbox()

from aircraftProperties import AircraftProperties
from wingboxGeometry import WingboxGeometry
import AerodynamicLoading

print(AircraftProperties.Fuselage["Fuselage length"])

wingbox = WingboxGeometry(forwardSpar=0.15, aftSpar=0.6)
print(wingbox.edgeCoordinates())
print(wingbox.calculateEnclArea())
wingbox.drawWingbox()

test = WingboxGeometry(forwardSpar=0.087654, aftSpar=0.87)
test.drawWingbox()

print(AerodynamicLoading.Dragacc(10))
AerodynamicLoading.drawgraphs()
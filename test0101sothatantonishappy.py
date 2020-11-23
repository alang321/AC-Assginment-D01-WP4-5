from math import sqrt, pi, radians, tan
from aircraftProperties import AircraftProperties

A = AircraftProperties.Planform["aspect ratio"]
Beta = sqrt(1-(AircraftProperties.Cruise_constants["mach at cruise"])**2)
eff = 0.95
sweep = AircraftProperties.Planform["half-chord sweep"]

y = (2*pi*A)/(2+sqrt(4+(((A*Beta)**2)/eff)*(1+((tan(radians(sweep))**2)/(Beta**2)))))
print(y)
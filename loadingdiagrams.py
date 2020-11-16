from aircraftProperties import AircraftProperties
import math

###---constants---###

#altitudes [m]
h_0 = 0                                                              #sea level 
rho_0 = 1.225                                                  
h_1 =     1                                                          # ???? (intermediate altitude)
rho_1 =  1
h_c = AircraftProperties.Cruise_constants["cruise altitude"]         #altitude at cruise 
rho_c = AircraftProperties.Cruise_constants["density at cruise"] 

#weights [N]
MTOW = AircraftProperties.Weight["MTOW"]                             #maximum takeoff weight
W_pl = AircraftProperties.Weight["payload at harmonic profile"]      #weight payload 
OEW = AircraftProperties.Weight["OEW"]                               #operating empty weight
ZFW = OEW + W_pl                                                     #zero fuel weight

#velocities
V_c = AircraftProperties.Cruise_constants["cruise velocity"]

#aircraft properties
S = AircraftProperties.Planform["surface area"]
C_L_max_flaps = AircraftProperties.Lift["CL max with flaps"]
C_L_max = AircraftProperties.Lift["CL max without flaps"]


###---velocity calculations---###

def get_V_S0(W, rho_h):                                                 #stall speed with flaps
    
    V_S0 = math.sqrt( (2 * W) / (rho_h * C_L_max_flaps * S) )
    
    return V_S0

def get_V_S1(W, rho_h):                                                 #stall speed without flaps
    
    V_S1 = math.sqrt( (2 * W) / (rho_h * C_L_max * S) )
    
    return V_S1






a = get_V_S0(MTOW, rho_0)

print(a)

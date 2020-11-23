from aircraftProperties import AircraftProperties
import math
import numpy as np
import matplotlib.pyplot as plt

###---constants---###

#altitudes [m]
h_0 = 0                                                              #sea level 
rho_0 = 1.225
T_0 = 288.15
h_1 =  4572                                                          # ???? (intermediate altitude)
rho_1 = 0.770816  
T_1 = 258.432
h_c = AircraftProperties.Cruise_constants["cruise altitude"]         #altitude at cruise 
rho_c = AircraftProperties.Cruise_constants["density at cruise"]
T_c = 226.733
g = 9.81
pi = 3.14159

#weights [N]
MTOW = AircraftProperties.Weight["MTOW"]                             #maximum takeoff weight
W_pl = AircraftProperties.Weight["payload at harmonic profile"]      #weight payload 
OEW = AircraftProperties.Weight["OEW"]                               #operating empty weight
ZFW = OEW + W_pl                                                     #zero fuel weight
weights_dic = {'OEW': OEW, 'ZFW': ZFW, 'MTOW': MTOW}

#velocities
M_c = AircraftProperties.Cruise_constants["mach at cruise"]


#aircraft properties
S = AircraftProperties.Planform["surface area"]
C_L_max_flaps = AircraftProperties.Lift["CL max with flaps"]
C_L_max = AircraftProperties.Lift["CL max without flaps"]
C_MAC = AircraftProperties.Planform["MAC"]


#maximimum load factors
n_max = 2.5
n_min = -1.0

#dictionaries
rho_h_dic = {'SL': rho_0, 'FL150': rho_1, 'FL310': rho_c}
T_dic = {'SL': T_0, 'FL150': T_1, 'FL310': T_c}
h_dic = {'SL': 0, 'FL150': 4572, 'FL310': 9449}
weights_dic = {'OEW': OEW, 'ZFW': ZFW, 'MTOW': MTOW}
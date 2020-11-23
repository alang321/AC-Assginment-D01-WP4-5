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
C_L_alpha = 4.62


#maximimum load factors
n_max = 2.5
n_min = -1.0

#dictionaries
rho_h_dic = {'SL': rho_0, 'FL150': rho_1, 'FL310': rho_c}
T_dic = {'SL': T_0, 'FL150': T_1, 'FL310': T_c}
h_dic = {'SL': 0, 'FL150': 4572, 'FL310': 9449}
weights_dic = {'OEW': OEW, 'ZFW': ZFW, 'MTOW': MTOW}

#quick calculations
Z_mo = h_c
R_1 = 1
R_2 = ZFW / MTOW
F_gz = 1 - (Z_mo/76200)
F_gm = math.sqrt(R_2 * math.tan(pi*R_1/4))
F_g = 0.5*(F_gz+F_gm)

###--- Gust Calculations ---###
def get_weight_factor(W, rho_h):

    weight_factor = (2 * W) / (S * rho_h * g * C_MAC * C_L_alpha)

    return weight_factor

def get_K_g(weight_factor):

    K_g = (0.88 * weight_factor) / (5.3 + weight_factor)

    return K_g

def get_U_ref(altitude):
    if altitude <= 4572:
        U_ref = 17.07 -0.000801*altitude
    else:
        U_ref = 17.07 - (0.000801*4572) - (0.000513998*(altitude-4572))
    return U_ref

def get_V_B(W, K_g, V_S1, V_C, U_ref):

    V_B = V_S1 * math.sqrt(1+ (K_g * rho_0 * U_ref * V_C * C_L_alpha)/(2*W))

    return V_B

def get_U_ds(U_ref, H):

    U_ds = U_ref * F_g * ((H/107)**(1/6))

    return U_ds

def get_U(U_ds, H, s):

    U = (U_ds/2)*(1-math.cos(pi*s/H))

    return U

def get_omega(V, H):

    omega = pi * V / H

    return omega

def time_constant(W, C_L_alpha, V, rho_h):

    time_constant = (W*2)/(S * C_L_alpha * V * rho_h * g)

    return time_constant

def gust_load_factor(t, U_ds, omega, time_constant):

    gust_load_factor = (U_ds/(2*g))*((omega*math.sin(omega*t))+((1/(1+(omega*time_constant**(-2))))*((math.exp(-t/time_constant))/time_constant)-(math.cos(omega*t)/time_constant)-(omega*math.sin(omega*t))))

    return gust_load_factor


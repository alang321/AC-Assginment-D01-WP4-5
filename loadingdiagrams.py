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




###---velocity calculations---###

def get_V_C(T_h):

    V_C = M_c * math.sqrt(T_h * 402.058)

    return V_C

def get_V_S0(W, rho_h):                                                 #stall speed with flaps
    
    V_S0 = math.sqrt( (2 * W) / (rho_h * C_L_max_flaps * S) )
    
    return V_S0

def get_V_S1(W, rho_h):                                                 #stall speed without flaps
    
    V_S1 = math.sqrt( (2 * W) / (rho_h * C_L_max * S) )
    
    return V_S1

def get_V_A(V_S1):

    V_A = math.sqrt(n_max)*V_S1

    return V_A

def get_V_D(V_C):

    V_D = V_C * 1.25

    return V_D

def get_V_F(rho_h):
    V_F1 = 1.6* get_V_S1(MTOW, rho_h) 
    V_F2 = 1.8* get_V_S1(ZFW, rho_h)
    V_F3 = 1.8* get_V_S0(ZFW, rho_h)

    if V_F1 > V_F2 and V_F1 > V_F3:
        return V_F1
    elif V_F2 > V_F3:
        return V_F2
    else:
        return V_F3
    


###---aerodynamic calculations---###
C_L_alpha = 4.62


###---gust calculations---###
Z_mo = h_c
R_1 = 1
R_2 = ZFW / MTOW
F_gz = 1 - (Z_mo/76200)
F_gm = math.sqrt(R_2 * math.tan(pi*R_1/4))
F_g = 0.5*(F_gz+F_gm)


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

def get_V_B(W, K_g, V_S1):

    V_B = V_S1 * math.sqrt(1+ (K_g * rho_0 * U_ref * V_C * C_L_alpha)/(2*W))

    return V_B

def get_U_ds(U_ref, H):

    U_ds = Uref * Fg * ((H/107)**(1/6))

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

    gust_load_factor = (U_ds/(2*g))*((omega*math.sin(omega*t))+((1/(1+(omega*time_constant**(-2))))*((e**(-t/time_constant))/time_constant)-(math.cos(omega*t)/time_constant)-(omega*math.sin(omega*t))))

    return gust_load_factor




###---speeds results---###

V_C_values = {}
V_S0_values = {}
V_S1_values = {}
V_D_values = {}
V_A_values = {}
V_B_values = {}
V_F_values = {}


          
for weights, values in weights_dic.items():
    for altitudes, temperatures in T_dic.items():

        nametagC = 'V_C {} {}'.format(weights, altitudes)
        nametagD = 'V_D {} {}'.format(weights, altitudes)

        V_C = get_V_C(temperatures)
        V_D = get_V_D(V_C)

        V_C_values[nametagC] = V_C
        V_D_values[nametagD] = V_D

    
    for altitudes, densitys in rho_h_dic.items():

        nametagS0 = 'V_S0 {} {}'.format(weights, altitudes)
        nametagS1 = 'V_S1 {} {}'.format(weights, altitudes)
        nametagA = 'V_A {} {}'.format(weights, altitudes)
        nametagB = 'V_B {} {}'.format(weights, altitudes)
        nametagF = 'V_F {} {}'.format(weights, altitudes)

        V_S0 = get_V_S0(values, densitys)
        V_S1 = get_V_S1(values, densitys)
        V_A = get_V_A(V_S1)
        V_F = get_V_F(densitys)

        altitude = h_dic[altitudes]
        U_ref = get_U_ref(altitude)
        weight_factor = get_weight_factor(values, densitys)
        K_g = get_K_g(weight_factor)
        V_B = get_V_B(values, K_g, V_S1)
        
        V_S0_values[nametagS0] = V_S0
        V_S1_values[nametagS1] = V_S1
        V_A_values[nametagA] = V_A
        V_B_values[nametagB] = V_B
        V_F_values[nametagF] = V_F
        


###---converting results to usable list---###        
        
V_all_dict = {}           
V_all_dict.update(V_A_values)
V_all_dict.update(V_C_values)
V_all_dict.update(V_D_values)
V_all_dict.update(V_F_values)
V_all_dict.update(V_S0_values)
V_all_dict.update(V_S1_values)
V_all_list = list(V_all_dict.values())



OEW_SL = []
OEW_FL150 = []
OEW_FL310 = []
ZFW_SL = []
ZFW_FL150 = []
ZFW_FL310 = []
MTOW_SL = []
MTOW_FL150 = []
MTOW_FL310 = []

for i in range(len(V_all_list)): 
    if i % 9 == 0:
        OEW_SL.append(V_all_list[i])
    elif i % 9 == 1:
        OEW_FL150.append(V_all_list[i])
    elif i % 9 == 2:
        OEW_FL310.append(V_all_list[i])
    elif i % 9 == 3:
        ZFW_SL.append(V_all_list[i])
    elif i % 9 == 4:
        ZFW_FL150.append(V_all_list[i])
    elif i % 9 == 5:
        ZFW_FL310.append(V_all_list[i])
    elif i % 9 == 6:
        MTOW_SL.append(V_all_list[i])
    elif i % 9 == 7:
        MTOW_FL150.append(V_all_list[i])
    elif i % 9 == 8:
        MTOW_FL310.append(V_all_list[i])

V_all_list_sorted = [OEW_SL, OEW_FL150, OEW_FL310, ZFW_SL, ZFW_FL150, ZFW_FL310, MTOW_SL, MTOW_FL150, MTOW_FL310]




###---MANEUVRE LOAD DIAGRAM---###

def f(x):
    return (x / V_S1)**2

def plot_maneuver(V_A, V_D, V_F, V_S0, V_S1):
                  
    speeds = [V_A, V_D, V_D, V_F, V_S1]  
    n_values = [n_max, n_max, 0, n_min, n_min]

    plt.plot(speeds, n_values, 'black')
    plt.title('Maneuvre Envelope')
    plt.xlabel('Velocity')
    plt.ylabel('Load')
    plt.text(V_A, n_max, 'V_A') #add text to diagram


    x1 = np.linspace(0,V_A, 1000)
    x2 = np.linspace(0, V_S1 * math.sqrt(2), 1000)
    x3 = np.linspace(0, V_S1, 1000)
    y2 = []      #flaps down curve n values

    for i in x2:
        a = (i / V_S0)**2
        a = min(a, 2)
        y2.append(a)

    plt.plot(x1, f(x1), 'black')  # (0,0) to V_A curve
    plt.plot(x2, y2, 'black')  #flaps down curve
    plt.plot(x3, -f(x3), 'black')
 
    return plt.show()

for sets in V_all_list_sorted:

    V_A = sets[0]
    V_D = sets[2]
    V_F = sets[3]
    V_S0 = sets[4]
    V_S1 = sets[5]

    plot_maneuver(V_A, V_D, V_F, V_S0, V_S1)




    
               
               


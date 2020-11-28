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
C_L_alpha = 4.0
AR = AircraftProperties.Planform["aspect ratio"]
half_c_sweep = AircraftProperties.Planform["half-chord sweep"]
eta = 0.95


#maximimum load factors
n_max = 2.5
n_min = -1.0

#dictionaries
altitudes = ['SL','FL150','FL310']
rho_h_dic = {'SL': rho_0, 'FL150': rho_1, 'FL310': rho_c}
T_dic = {'SL': T_0,'FL150': T_1,'FL310': T_c}
h_dic = {'SL': 0,'FL150': 4572,'FL310': 9449}
weights_dic = {'OEW': OEW, 'ZFW': ZFW, 'MTOW': MTOW}

#quick calculations
Z_mo = h_c
R_1 = 1
R_2 = ZFW / MTOW
F_gz = 1 - (Z_mo/76200)
F_gm = math.sqrt(R_2 * math.tan(pi*R_1/4))
F_g_init = 0.5*(F_gz+F_gm)

###---maneuver calculations---###

def get_VC_VD(T_h):

    V_C = M_c * math.sqrt(T_h * 402.058)
    V_D = V_C * 1.25

    return V_C, V_D

def get_V_S0(W, rho_h):                                                 #stall speed with flaps
    
    V_S0 = math.sqrt( (2 * W) / (rho_h * C_L_max_flaps * S) )
    
    return V_S0

def get_V_S1(W, rho_h):                                                 #stall speed without flaps
    
    V_S1 = math.sqrt( (2 * W) / (rho_h * C_L_max * S) )
    
    return V_S1

def get_V_A(V_S1):

    V_A = math.sqrt(n_max)*V_S1

    return V_A

def get_V_F(rho_h):
    V_F1 = 1.6* get_V_S1(MTOW, rho_h) 
    V_F2 = 1.8* get_V_S1(ZFW, rho_h)
    V_F3 = 1.8* get_V_S0(ZFW, rho_h)

    if V_F1 > V_F2 and V_F1 > V_F3:
        V_F = V_F1
    elif V_F2 > V_F3:
        V_F = V_F2
    else:
        V_F = V_F3

    return V_F

###--- Gust Calculations ---###
def get_C_L_alpha_h(V, T_h):

    M = V/(math.sqrt(T_h * 402.058))

    if M < 0.70:
        beta = math.sqrt(1-M**2)
    else:
        beta = math.sqrt(1-0.70**2)

    term1 = 2*pi*AR
    term2 = 4+(AR*beta/eta)**2
    term3 = ((math.tan(half_c_sweep))/beta)**2


    C_L_alpha_h = term1/(2 + math.sqrt(4+(term2)*(term3)))

    return C_L_alpha_h


def get_F_g(altitude):

    F_g = ((1.0 - F_g_init)/Z_mo) * altitude + F_g_init

    return F_g

def get_weight_factor(W, rho_h, C_L_alpha_h):

    weight_factor = (2 * W) / (S * rho_h * g * C_MAC * C_L_alpha_h)

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

def get_V_B(W, K_g, V_S1, V_C, U_ref, C_L_alpha_h):

    V_B = V_S1 * math.sqrt(1+ (K_g * rho_0 * U_ref * V_C * C_L_alpha_h)/(2*W))

    return V_B

def get_U_ds(U_ref, F_g, H):

    U_ds = U_ref * F_g * ((H/107)**(1/6))

    return U_ds

def get_U(U_ds, H, s):

    U = (U_ds/2)*(1-math.cos(pi*s/H))

    return U

def get_omega(V, H):

    omega = pi * V / H

    return omega

def get_time_constant(W, V, rho_h, C_L_alpha_h):

    time_constant = (W*2)/(S * C_L_alpha_h * V * rho_h * g)

    return time_constant

def get_gust_load_factor(t, U_ds, omega, time_constant):
    term1 = U_ds/(2*g)
    term2 = omega * math.sin(omega*t)
    term3 = 1 / (1+(omega*time_constant)**(-2))
    term4 = (1/time_constant)* math.exp(-t/time_constant)
    term5 = math.cos(omega*t)/time_constant

    gust_load_factor = term1 * (term2 + term3*(term4-term5-term2))

    return gust_load_factor


###---Find Maximum Load Factor Conditions ---###
max_load_factor = 0
V_max_load = 0
W_max_load = 0
H_max_load = 0
h_max_load = 0
t_max_load = 0
U_ds_max_load = 0
omega_max_load = 0
n_max_load = 0
time_constant_max_load = 0

for weights, values in weights_dic.items():
    for keys in altitudes:
        V_tot = []

        h = h_dic[keys]
        T_h = T_dic[keys]
        rho_h = rho_h_dic[keys]
    
        U_ref = get_U_ref(h)
        F_g = get_F_g(h)
        weight_factor = get_weight_factor(values, rho_h, C_L_alpha)
        K_g = get_K_g(weight_factor)
        
        V_C, V_D = get_VC_VD(T_h)
        V_S0 = get_V_S0(values, rho_h)
        V_S1 = get_V_S1(values, rho_h)
        V_A = get_V_A(V_S1)
        V_F = get_V_F(rho_h)
        V_B = get_V_B(values, K_g, V_S1, V_C, U_ref, C_L_alpha)
        
        V_tot.append(V_C)
        V_tot.append(V_D)
        V_tot.append(V_S0)
        V_tot.append(V_S1)
        V_tot.append(V_A)
        V_tot.append(V_F)
        V_tot.append(V_B)

        for V in V_tot:
            C_L_alpha_h = get_C_L_alpha_h(V, T_h)
            time_constant = get_time_constant(values, V, rho_h, C_L_alpha_h)
            
            for H in range(9, 108, 1):
                omega = get_omega(V, H)
                U_ds0 = get_U_ds(U_ref, F_g, H)
                n = 2*pi / omega

                if V is not V_D:
                    U_ds = U_ds0
                else:
                    U_ds = U_ds0 * 0.5

                for t in np.linspace(0, n, 100):
                    
                    gust_load_factor = get_gust_load_factor(t, U_ds, omega, time_constant)

                    if gust_load_factor > max_load_factor:                      
                        max_load_factor = gust_load_factor
                        V_max_load = V
                        W_max_load = weights
                        H_max_load = H
                        h_max_load = keys
                        t_max_load = t
                        n_max_load = n

                        U_ds_max_load = U_ds
                        omega_max_load = omega
                        time_constant_max_load = time_constant
                        
                        
print('Maximum gust load factor: {}'.format(max_load_factor))
print('Occurst at:')
print('Weight: {}'.format(W_max_load))
print('Altitude: {}'.format(h_max_load))
print('Velocity: {}'.format(V_max_load))
print('Gust gradient distance: {}'.format(H_max_load))

print(n_max_load)
###--- Gust Loads Plotters ---###
                              
def gust_load_plotter(n, U_ds, omega, time_constant):
    def f(x):
        return get_gust_load_factor(x, U_ds, omega, time_constant)

    x1 = np.linspace(0, n, 100)
    y1 = []

    plt.xlabel('t')
    plt.ylabel('n')
    plt.title('Plot of the variation of the gust load factor as function of time')

    for x in x1:
        y1.append(f(x))

    plt.plot(x1, y1, 'black')

    return plt.show()

gust_load_plotter(n_max_load, U_ds_max_load, omega_max_load, time_constant_max_load)

###--- V/n diagram ---###

def get_V_n_diagram(V, weight, altitude):
    value = weights_dic[weight]
    max_load_factor = 0
    h = h_dic[altitude]
    T_h = T_dic[altitude]
    rho_h = rho_h_dic[altitude]
        
    U_ref = get_U_ref(h)
    F_g = get_F_g(h)
    V_C, V_D = get_VC_VD(T_h)

    C_L_alpha_h = get_C_L_alpha_h(V, T_h)
    time_constant = get_time_constant(value, V, rho_h, C_L_alpha_h)
            
    for H in range(9, 108, 3):
        omega = get_omega(V, H)
        U_ds0 = get_U_ds(U_ref, F_g, H)
        n = 2*pi / omega

        if V <= V_C:
            U_ds = U_ds0
        else:
            U_ds = U_ds0 * 0.5

        for t in np.linspace(0, n, 100):
                        
            gust_load_factor = get_gust_load_factor(t, U_ds, omega, time_constant)

            if gust_load_factor > max_load_factor:                      
                max_load_factor = gust_load_factor

    return max_load_factor

print(get_VC_VD(T_dic['SL']))
print(get_V_n_diagram(327, 'OEW', 'SL'))



def plot_V_n_diagram(weight, altitude):
    T_h = T_dic[altitude]
    V_C, V_D = get_VC_VD(T_h)

    x1 = np.linspace(0.01, V_C, 100)
    x2 = np.linspace(V_C, V_D, 50)
    

    y1 = []
    y2 = []
    y3 = np.linspace(1+get_V_n_diagram(V_C, weight, altitude), 1+get_V_n_diagram(V_D, weight, altitude), 50)
    y4 = np.linspace(1 - get_V_n_diagram(V_C, weight, altitude), 1 - get_V_n_diagram(V_D, weight, altitude), 50)
    
    plt.xlabel('V')
    plt.ylabel('n')
    plt.title('Plot of the variation of the gust load factor as function of time')

    for x in x1:
        y1.append(get_V_n_diagram(x, weight, altitude)+1)
        y2.append(1-get_V_n_diagram(x, weight, altitude))
        

    plt.plot(x1, y1, 'black')
    plt.plot(x1, y2, 'black')
    plt.plot(x2, y3, 'black')
    plt.plot(x2, y4, 'black')
    plt.vlines(V_D, (1-get_V_n_diagram(V_D, weight, altitude)),(get_V_n_diagram(V_D, weight, altitude)+1),  color = 'black',)

    return plt.show()

plot_V_n_diagram('OEW', 'SL')
print(1-get_V_n_diagram(V_D, 'OEW', 'SL'))


                        

                        

        
            
        


        
       






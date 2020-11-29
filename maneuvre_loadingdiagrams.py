from aircraftProperties import AircraftProperties
import math
import numpy as np
import matplotlib.pyplot as plt
import xlsxwriter 

###---Create xlsx file---###

workbook = xlsxwriter.Workbook('Load_Table.xlsx') 
worksheet = workbook.add_worksheet() 
  
#create columns
worksheet.write('A1', 'Case') 
worksheet.write('B1', 'Design Speed') 
worksheet.write('C1', 'Speed Value') 
worksheet.write('D1', 'Weight Key') 
worksheet.write('E1', 'Weight Value')
worksheet.write('F1', 'Altitude Key')
worksheet.write('G1', 'Altitude [m]')
worksheet.write('H1', 'n') 
worksheet.write('I1', 'Fuel Weight') 
worksheet.write('J1', 'Thrust') 

global worksheet_key
worksheet_key = 2    



###---constants---###

#altitudes [m]
h_0 = 0                                                              #sea level 
rho_0 = 1.225
T_0 = 288.15
h_1 =  6096                                                         # ???? (intermediate altitude)
rho_1 = 0.652694 
T_1 = 248.526
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
W_fuel = MTOW - ZFW

weights_dic = {'OEW': OEW, 'ZFW': ZFW, 'MTOW': MTOW}
TOW0_5 = OEW + 0.5*W_pl + W_fuel
TOW0 = OEW + W_fuel

Cff = (0.99*0.99*0.995*0.98*0.856)
CW1 = Cff*W_fuel + OEW + W_pl
CW0_5 = Cff*W_fuel + OEW + 0.5*W_pl
CW0 = Cff*W_fuel + OEW

Loiff = Cff*0.99*0.988
LoiW1 = Loiff*W_fuel + OEW + W_pl
LoiW0_5 = Loiff*W_fuel + OEW + 0.5*W_pl
LoiW0 = Loiff*W_fuel + OEW

Lff = Loiff*0.965*0.992
LW1 = Lff*W_fuel + OEW + W_pl
LW0_5 = Lff*W_fuel + OEW + 0.5*W_pl
LW0 = Lff*W_fuel + OEW

weights_altitudes = {'MTOW': ['SL', MTOW],'TOW_0.5': ['SL', TOW0_5], 'TOW_0':['SL', TOW0]
                     , 'CW_1': ['FL310', CW1], 'CW_0.5': ['FL310', CW0_5], 'CW_0': ['FL310', CW0],
                     'LoiW_1': ['FL200', LoiW1], 'LoiW_0.5': ['FL200', LoiW0_5], 'LoiW_0': ['FL200', LoiW0],
                     'LW_1': ['SL', LW1], 'LW_0.5': ['SL', LW0_5], 'LW_0': ['SL', LW0],}


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
rho_h_dic = {'SL': rho_0, 'FL200': rho_1, 'FL310': rho_c}
T_dic = {'SL': T_0, 'FL200': T_1, 'FL310': T_c}
h_dic = {'SL': 0, 'FL200': 4572, 'FL310': 9449}
weights_dic = {'OEW': OEW, 'ZFW': ZFW, 'MTOW': MTOW}

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

          
for weights, values in weights_altitudes.items():

    nametagC = 'V_C {} {}'.format(weights, values[0])
    nametagD = 'V_D {} {}'.format(weights, values[0])
    nametagS0 = 'V_S0 {} {}'.format(weights, values[0])
    nametagS1 = 'V_S1 {} {}'.format(weights, values[0])
    nametagA = 'V_A {} {}'.format(weights, values[0])
    nametagB = 'V_B {} {}'.format(weights, values[0])
    nametagF = 'V_F {} {}'.format(weights, values[0])

    T_h = T_dic[values[0]]
    densitys = rho_h_dic[values[0]]
        
    V_C, V_D = get_VC_VD(T_h)
    V_S0 = get_V_S0(values[1], densitys)
    V_S1 = get_V_S1(values[1], densitys)
    V_A = get_V_A(V_S1)
    V_F = get_V_F(densitys)

    altitude = h_dic[values[0]]
    U_ref = get_U_ref(altitude)
    weight_factor = get_weight_factor(values[1], densitys)
    K_g = get_K_g(weight_factor)
    V_B = get_V_B(values[1], K_g, V_S1)

    V_C_values[nametagC] = V_C
    V_D_values[nametagD] = V_D 
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



MTOW = []
TOW_0_5 = []
TOW_0 = []
CW_1 = []
CW_0_5 = []
CW_0 = []
LoiW_1 = []
LoiW_0_5 = []
LoiW_0 = []
LW_1 = []
LW_0_5 = []
LW_0 = []

for i in range(len(V_all_list)): 
    if i % 12 == 0:
        MTOW.append(V_all_list[i])
    elif i % 12 == 1:
        TOW_0_5.append(V_all_list[i])
    elif i % 12 == 2:
        TOW_0.append(V_all_list[i])
    elif i % 12 == 3:
        CW_1.append(V_all_list[i])
    elif i % 12 == 4:
        CW_0_5.append(V_all_list[i])
    elif i % 12 == 5:
        CW_0.append(V_all_list[i])
    elif i % 12 == 6:
        LoiW_1.append(V_all_list[i])
    elif i % 12 == 7:
        LoiW_0_5.append(V_all_list[i])
    elif i % 12 == 8:
        LoiW_0.append(V_all_list[i])
    elif i % 12 == 9:
        LW_1.append(V_all_list[i])
    elif i % 12 == 10:
        LW_0_5.append(V_all_list[i])
    elif i % 12 == 11:
        LW_0.append(V_all_list[i])

V_all_list_sorted = [MTOW, TOW_0_5, TOW_0, CW_1, CW_0_5, CW_0, LoiW_1, LoiW_0_5, LoiW_0, LW_1, LW_0_5, LW_0]
print(len(V_all_list_sorted))

print(V_all_list_sorted[0])

print(LW0)
print(LW0_5)
print(LW1)
print(OEW)
print(ZFW)



###---MANEUVRE LOAD DIAGRAMS---###

def plot_maneuver(V_A, V_D, V_F, V_S0, V_S1, title = 'Manoeuvre diagram'):

    plt.figure()
    
    def f(x):
        return (x / V_S1)**2
                      
    speeds = [V_A, V_D, V_D, V_F, V_S1]  
    n_values = [n_max, n_max, 0, n_min, n_min]

    plt.plot(speeds, n_values, 'black')
    plt.title(title)
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Load (n)')
    
    
    plt.axvline(x=V_S0, color = 'purple', linestyle = '--', label = "V_S0")
    
    
    #plt.text(V_S1+5, -1.25, 'V_S1') #add text to diagram 
    plt.axvline(x=V_S1, color = 'red', linestyle = '--', label = "V_S1")

    #plt.text(V_A+5, -1.25, 'V_A') #add text to diagram 
    plt.axvline(x=V_A, color = 'green', linestyle = '--', label = "V_A")
    

    #plt.text(V_F+5, -1.25, 'V_C') #add text to diagram 
    plt.axvline(x=V_F, color = 'orange', linestyle = '--', label = "V_F")
    
    plt.axvline(x=V_C, color = 'aqua', linestyle = '--', label = "V_C")

    
    #plt.text(V_D+5, -1.25, 'V_D') #add text to diagram 
    plt.axvline(x=V_D, color = 'navy', linestyle = '--', label = "V_D")

   
    

    
    plt.text(V_S1 * math.sqrt(2) - 70, 2.1, 'Flaps') #add text to diagram 
    plt.text(V_S1 * math.sqrt(2) - 70, 1.8, 'Down') #add text to diagram 


    plt.axhline(y=0, color = 'black', linestyle = ':')
    
    plt.xticks(np.arange(0, 450, 50))  #fixing the x axis 
    plt.ylim(-1.5, 3)
    plt.xlim(0,450)
    
    plt.yticks(np.arange(-1.5, 3, 0.5))  #fixing the y axis 


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
    
    leg = plt.legend();

    plt.show()
    
    ###---updating excel file---###
    
    global worksheet_key
    a = int(( worksheet_key - 2 ) / 5)
    
    speeds_dic = {"V_A": V_A, "V_D": V_D, "V_F": V_F, "V_S0": V_S0, "V_S1": V_S1}
    
    fuel_weights = {'MTOW': W_fuel,'TOW_0.5': W_fuel, 'TOW_0':W_fuel, 'CW_1': Cff*W_fuel, 
                    'CW_0.5': Cff*W_fuel, 'CW_0': Cff*W_fuel,'LoiW_1': Loiff*W_fuel, 
                    'LoiW_0.5': Loiff*W_fuel, 'LoiW_0': Loiff*W_fuel,'LW_1': Lff*W_fuel, 
                    'LW_0.5': Lff*W_fuel, 'LW_0': Lff*W_fuel,}
    
    loads = {"V_A": n_max, "V_D": n_max, "V_F": n_min, "V_S0": ( V_S0 / V_S1 ) ** 2, "V_S1": 1}

    for i in range(worksheet_key, worksheet_key+5):
        
        worksheet.write("A"+str(i),str(i-1))
        worksheet.write("B"+str(i),str(list(speeds_dic.keys())[i-worksheet_key]))
        worksheet.write("C"+str(i),str(list(speeds_dic.values())[i-worksheet_key]))
        worksheet.write("D"+str(i),list(weights_altitudes.keys())[a])
        worksheet.write("E"+str(i),str(list(weights_altitudes.values())[a][1]))
        worksheet.write("F"+str(i),str(list(weights_altitudes.values())[a][0]))
        worksheet.write("G"+str(i),h_dic[str(list(weights_altitudes.values())[a][0])])
        worksheet.write("I"+str(i),str(fuel_weights[str(list(weights_altitudes.keys())[a])]))
        worksheet.write("H"+str(i),str(loads[list(speeds_dic.keys())[i-worksheet_key]]))
        
        
    worksheet_key += 5
    

for i in range(len(V_all_list_sorted)):

    V_A = V_all_list_sorted[i][0]
    V_D = V_all_list_sorted[i][2]
    V_F = V_all_list_sorted[i][3]
    V_S0 = V_all_list_sorted[i][4]
    V_S1 = V_all_list_sorted[i][5]
        
    if i  == 0:
        plot_maneuver(V_A, V_D, V_F, V_S0, V_S1, title = 'Manoeuvre envelope at MTOW and SL')
    elif i % 12 == 1:
        plot_maneuver(V_A, V_D, V_F, V_S0, V_S1, title = 'Manoeuvre envelope at TOW with half payload and SL')
    elif i % 12 == 2:
        plot_maneuver(V_A, V_D, V_F, V_S0, V_S1, title = 'Manoeuvre envelope at TOW with zero payload and SL')
    elif i % 12 == 3:
        plot_maneuver(V_A, V_D, V_F, V_S0, V_S1, title = 'Manoeuvre envelope at CW with full payload and FL310')
    elif i % 12 == 4:
        plot_maneuver(V_A, V_D, V_F, V_S0, V_S1, title = 'Manoeuvre envelope at CW with half payload and FL310')
    elif i % 12 == 5:
        plot_maneuver(V_A, V_D, V_F, V_S0, V_S1, title = 'Manoeuvre envelope at CW with zero payload and FL310')
    elif i % 12 == 6:
        plot_maneuver(V_A, V_D, V_F, V_S0, V_S1, title = 'Manoeuvre envelope at LoiW with full payload and FL310')
    elif i % 12 == 7:
        plot_maneuver(V_A, V_D, V_F, V_S0, V_S1, title = 'Manoeuvre envelope at LoiW with half payload and FL310')
    elif i % 12 == 8:
        plot_maneuver(V_A, V_D, V_F, V_S0, V_S1, title = 'Manoeuvre envelope at LoiW with zero payload and FL310')
    elif i % 12 == 9:
        plot_maneuver(V_A, V_D, V_F, V_S0, V_S1, title = 'Manoeuvre envelope at LW with full payload and FL310')
    elif i % 12 == 10:
        plot_maneuver(V_A, V_D, V_F, V_S0, V_S1, title = 'Manoeuvre envelope at LW with half payload and FL310')
    elif i % 12 == 11:
        plot_maneuver(V_A, V_D, V_F, V_S0, V_S1, title = 'Manoeuvre envelope at LW with zero payload and FL310')


workbook.close() 

     
               


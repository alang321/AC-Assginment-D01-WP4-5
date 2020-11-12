from aircraftProperties import AircraftProperties

###---constants---###

#altitudes [m]
h_0 = 0         #sea level 
h_1 =   2       # ????
h_c = 9448      #altitude at cruise 

#weights [N]
W_pl = 66000 * 9.81                         #weight payload 
OEW = AircraftProperties.Weight["OEW"]      #operating empty weight
ZFW = OEW + W_pl                            #zero fuel weight
MTOW = AircraftProperties.Weight["MTOW"]    #maximum takeoff weight

#V_s0 =  
#V_s1 = 
#V_A = 
#V_C = 
#V_D = 

print(OEW)


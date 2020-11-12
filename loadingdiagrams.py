from aircraftProperties import AircraftProperties

###---constants---###

#altitudes [m]
h_0 = 0         #sea level 
h_1 =   2       # ????
h_c = AircraftProperties.Cruise_constants["cruise altitude"]   #altitude at cruise 

#weights [N]
MTOW = AircraftProperties.Weight["MTOW"]    #maximum takeoff weight
W_pl = AircraftProperties.Weight["payload at harmonic profile"]                #weight payload 
OEW = AircraftProperties.Weight["OEW"]      #operating empty weight
ZFW = OEW + W_pl                            #zero fuel weight

#V_s0 =  
#V_s1 = 
#V_A = 
#V_C = 
#V_D = 

print(h_c)


g = 9.80665
R = 287.0
T0 = 288.15
p = 101325.0
a = -0.0065


ch=float(input('\nEnter your choice:  '))

if (ch == 1):

    print('\n    **** ISA calculator troposhere ****')

    
    h=float(input('\nEnter altitude [m]    '))
        
    if h <= 11000.0:
        h1=h
        
        T1= T0 + (a*h)
        C1=T1 - 273.15
        print('\nTemperature:',round(T1,2),'K', '(',round(C1,1),'C)')
    
        
        p1= p*((T1/T0)**(-(g/(a*R))))
        Sl= (p1/p)*100
        print('Pressure:',round(p1), '(',round(Sl),'% SL)')
    
        
        rho= p1/(R*T1)
        Slr= (rho/1.225)*100
        print('density:',round(rho,4),'(',round(Slr),'% SL)')
        print('\nReady.')
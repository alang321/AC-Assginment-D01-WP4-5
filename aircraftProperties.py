class AircraftProperties:
    Weight = {
        "MTOW": 2014954.46240041,  # $N$
        "OEW": 891618.429884555,  # $N$
        "fuel fraction": 0.23394240822931,  # -
        "Design TOW": 1866756.29625086,  # $N$
        "design fuel fraction": 0.309146329830603,  # $N$
        "primary components of OEW": 80.4532673513519,  # %
        "secondary components of OEW": 8.95422584993104,  # %
        "tertiary components of OEW": 10.5925067987171,  # %
        "payload at harmonic profile": 66000 * 9.81  # $N$
    }
    Planform = {
        "surface area": 362.728076039678,  # $m^{2}$
        "quarter-chord sweep": 28.6796,  # $deg$
        "half-chord sweep": 25.9829,  # $deg$
        "leading edge sweep": 31.4151,  # $deg$
        "taper ratio": 0.2996,  # -
        "aspect ratio": 9,  # -
        "span": 57.1362641792155,  # $m$
        "root chord": 9.77,  # $m$
        "tip chord": 2.927092,  # $m$
        "MAC": 6.96319316876988,  # $m$
        "xMAX": 7.04,  # $m$
        "yMAC": 11.72,  # $m$
        "dihedral": 2.123,  # $deg$
        "forward spar fraction of chord": 0.15,  # -
        "aft spar fraction of chord": 0.6,  # -
        "cd0 wing": 0.00317965138,
        "e": 0.625
    }
    Airfoil = {
        "t/c": 14.22,  # $\%$
        "airfoil identifier": "NASA20714",  # -
        "zero lift aoa": -4.75 #deg
    }
    Aileron = {
        "span fract of wingspan": 0.28,  # -
        "chord ratio": 0.3,  # -
        "max deflection angle": 35,  # $deg$
        "inboard edge fract of wing span": 0.7  # -
    }
    High_Lift_Devices = {
        "max lift coefficient": 2.4,  # -
        "inboard edge TEHLD": 3.244925,  # m
        "outboard edge TEHLD": 20.0445160995897,  # m
        "inboard edge LEHLD": 4.68427398607778,  # m
        "outboard edge LEHLD": 26.3442211594607,  # m
        "TEHLD chord ratio": 0.3,  # -
        "LEHLD chord ratio": 0.3,  # -
        "Required scrape angle": 11.8  # deg
    }
    Engine = {
        "thrust": 431.5,  # $kN$
        "bypass ratio": 9.6,  # -
        "weight": 74066,  # $N$
        "specific fuel consumption": 0.478  # $lbs/lbf/h$
    }
    Empennage = {
        "Horizontal tail volume": 1,  # $m^3$
        "Vertical tail volume": 0.09,  # $m^3$
        "moment arm horizontal tail": 34.1485952096708,  # m
        "moment arm vertical tail": 34.1485952096708,  # m
        "Surface area horizontal tail": 73.9633840189491,  # $m^2$
        "Aspect ratio horizontal tail": 4.5,  # -
        "sweep angle at quarter cord horizontal": 33.6796,  # deg
        "taperratio horizontal tail": 0.37,  # -
        "span horizontal tail": 18.2437723096204,  # m
        "horizontal tail cord root": 5.91849872169355,  # m
        "horizontal tail chord tip": 2.18984452702661,  # m
        "mean aerodynamic cord horizontal tail": 4.33994439110463,  # m
        "y position horizontal tail": 3.86182041590504,  # m
        "Surface area vertical tail": 54.6213814814763,  # $m^2$
        "Aspect vertical tail": 1.5,  # -
        "sweep angle at quarter cord vertical": 38.7,  # deg
        "taperratio vertical tail": 0.48,  # -
        "span vertical tail": 9.05163367697868,  # m
        "vertical tail cord root": 8.15462493421502,  # m
        "vertical tail cord tip": 3.91421996842321,  # m
        "mean aerodynamic cord vertical tail": 6.28273445382044,  # m
        "y position vertical tail": 1.99788310888268  # m
    }
    Fuselage = {
        "Fuselage length": 75.22376,  # m
        "Fuselage diameter": 6.48985,  # m
        "Tailcone": 21.416505,  # m
        "Nosecone": 8.1123125,  # m
        "Fuselage thickness": 0.179924999999999,  # m
        "passenger number max": 520,  # -
        "pressurised volume": 7182.24188198783,  # $m^3$
        "luggage volume": 69.7411764705882,  # $m^3$
        "seats abreast": 10  # -
    }
    Landing_gear = {
        "wheels main landing gear": 8,  # -
        "wheels nose landing gear": 2,  # -
        "nose gear length": 4.350002349,  # m
        "main gear length": 4.350002349  # m
    }
    CG = {
        "most forward center of gravity": 33.0952561849258,  # m
        "most aft center of gravity": 33.5527887903292  # m
    }
    Fuel = {
        "total fuel volume required": 88.898,  # $m^3$
        "fuel density": 804, # kg/m^3
    }
    Cruise_constants = {
        "cruise altitude": 9448,             # $m$
        "mach at cruise": 0.77,
        "Temp at cruise": 226.738,  # $C$
        "speed of sound at cruise": 301.833, # $ms^{-1}$
        "cruise velocity": 232.412, # $ms^{-1}$
        "density at cruise": 0.4436323254,
    }
    Lift = {
        "CL max with flaps": 2.62,
        "CL max without flaps": 1.631
    }
    WingboxMaterial = {
        "Aluminum Name": "AL6061-T",
        "shear modulus": 26 * 10**9,
        "e modulus": 68.9 * 10**9,
        "density": 2700, #kg/m^3
        "yield strength": 2.76 * 10**8, #Pa
        "edge crack strength": 3.27 * 10**8, #Pa
        "center crack strength": 2.92 * 10**8, # Pa
        "ultimate strength": 3.10 * 10**8, #Pa
        "shear strength": 2.07 * 10**8, #Pa
        "poisson's ratio": 0.33 
    }
    Safety = {
        "ExternalLoadSF": 1.5
    }
    Rivets = {
        "column minimum rivet pitch": 30, #mm, assumed value
        "c": 2.1, #mm
    }


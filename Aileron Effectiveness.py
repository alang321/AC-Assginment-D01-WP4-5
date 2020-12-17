import matplotlib.pyplot as plt

SurfaceArea = 362.73
chord = 3.07
ShearModulus = 2.60E10
K = 8935860.634              # value of wingbox from WP4
# dCL_dXi = 4.9562
# dCm_dXi = -0.7735
dCL_dAlpha = 6.039
e = 0.125

# Sea Level:
rho = 1.225
AileronEffSea = []
SpeedValuesSea = []
for speed in range(300):
    dCL_dXi = 2.38 + 0.0241 * speed - 2.81E-5 * (speed ** 2)
    dCm_dXi = -0.419 - 3.7E-3 * speed + 7.67E-6 * (speed ** 2)
    Top = 0.5 * rho * (speed**2) * SurfaceArea * chord * dCm_dXi * dCL_dAlpha + K * dCL_dXi
    Bot = (K - 0.5 * rho * (speed**2) * SurfaceArea * chord * e * dCL_dAlpha) * dCL_dXi
    SpeedValuesSea.append(speed)
    AileronEffSea.append(Top/Bot)
    if Top/Bot <= 0:
        break
# Cruise:
rho = 0.4416
AileronEffCruise = []
SpeedValuesCruise = []
for speed in range(300):
    dCL_dXi = 3.57 + 9.01E-3 * speed + 1.28E-5 * (speed ** 2)
    dCm_dXi = -0.782 + 1.22E-3 * speed - 6.39E-6 * (speed ** 2)
    Top = 0.5 * rho * (speed**2) * SurfaceArea * chord * dCm_dXi * dCL_dAlpha + K * dCL_dXi
    Bot = (K - 0.5 * rho * (speed**2) * SurfaceArea * chord * e * dCL_dAlpha) * dCL_dXi
    SpeedValuesCruise.append(speed)
    AileronEffCruise.append(Top/Bot)
    if Top/Bot <= 0:
        break

# graph plotting
plt.plot(SpeedValuesSea, AileronEffSea, label='Sea Level')
plt.plot(SpeedValuesCruise, AileronEffCruise, label='Cruise altitude')
plt.axvline(x=117.9, color = 'red', linestyle = '--', label = "V_maneuvre")
plt.axvline(x=191.67, color = 'red', linestyle = '--')
plt.title('Aileron Effectiveness at Sea Level and Cruise Altitude vs. Velocity')
plt.xlabel('Velocity [m/s]')
plt.ylabel('Aileron Effectiveness [-]')
plt.ylim(0, 1)
plt.xlim(0)
plt.grid(True)
plt.legend()
plt.show()

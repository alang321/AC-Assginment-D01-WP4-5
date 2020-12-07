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
    if speed <= 137.9:
        dCL_dXi = 4.956167994
        dCm_dXi = -0.7735059875
    elif speed <= 157.90:
        dCL_dXi = 4.984816364
        dCm_dXi = -0.7448576176
    elif speed <= 177.90:
        dCL_dXi = 4.956167994
        dCm_dXi = -0.7735059875
    elif speed <= 197.90:
        dCL_dXi = 4.956167994
        dCm_dXi = -0.7735059875
    elif speed <= 217.90:
        dCL_dXi = 4.956167994
        dCm_dXi = -0.7735059875
    elif speed <= 237.90:
        dCL_dXi = 4.956167994
        dCm_dXi = -0.7735059875
    elif speed <= 257.90:
        dCL_dXi = 4.956167994
        dCm_dXi = -0.7735059875
    elif speed <= 277.90:
        dCL_dXi = 4.956167994
        dCm_dXi = -0.7735059875
    elif speed <= 297.90:
        dCL_dXi = 4.956167994
        dCm_dXi = -0.7735059875
    else:
        dCL_dXi = 4.956167994
        dCm_dXi = -0.7735059875
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
    Top = 0.5 * rho * (speed**2) * SurfaceArea * chord * dCm_dXi * dCL_dAlpha + K * dCL_dXi
    Bot = (K - 0.5 * rho * (speed**2) * SurfaceArea * chord * e * dCL_dAlpha) * dCL_dXi
    SpeedValuesCruise.append(speed)
    AileronEffCruise.append(Top/Bot)
    if Top/Bot <= 0:
        break

# graph plotting
plt.plot(SpeedValuesSea, AileronEffSea, label='Sea Level')
plt.plot(SpeedValuesCruise, AileronEffCruise, label='Cruise')
plt.title('Aileron Effectiveness vs Velocity')
plt.xlabel('Velocity [m/s]')
plt.ylabel('Aileron Effectiveness [-]')
plt.ylim(0, 1)
plt.xlim(0)
plt.grid(True)
plt.legend()
plt.show()

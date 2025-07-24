import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ML8
pressure = np.array([10.18, 18.81, 24.84, 30.29, 34.34])
a_values = np.array([4.4702, 4.4304, 4.4174, 4.3963, 4.38])
a_errors = np.array([0.0011, 0.013, 0.0027, 0.0011, 0.0054])

# ML3
pressure1 = np.array([3.68, 9.85, 15.08, 20.21, 25.51, 31.06, 35.71])
a_values1 = np.array([4.4941, 4.4652, 4.4479, 4.4378, 4.4205, 4.4071, 4.3911])
a_errors1 = np.array([0.0061, 0.0031, 0.0035, 0.0068, 0.0046, 0.0045, 0.0075])

# Convert to volumes
a_values_cubed = a_values ** 3
a_errors_cubed = 3 * (a_values ** 2) * a_errors

a_values_cubed1 = a_values1 ** 3
a_errors_cubed1 = 3 * (a_values1 ** 2) * a_errors1

# ZrC theoretical parameters
Vo = 104.75
B = 217.7
Bo = 3.934

# Birch-Murnaghan equation
def birch_murnaghan(V, V0, B0, B0_prime):
    eta = (V0 / V)
    return ((3 / 2) * B0) * (eta**(7/2) - eta**(5/2)) * (1 + (3 / 4)*(B0_prime - 4)*(eta**(2/3) - 1))

# Curve fitting
initial_guess = [a_values_cubed[0], 100, 4]
popt1, _ = curve_fit(birch_murnaghan, a_values_cubed, pressure, p0=initial_guess)
popt2, _ = curve_fit(birch_murnaghan, a_values_cubed1, pressure1, p0=initial_guess)

# Define plotting ranges
volume_range_ZrC = np.linspace(Vo * 0.9, Vo * 1.1, 200)
pressure_ZrC_fit = birch_murnaghan(volume_range_ZrC, Vo, B, Bo)

volume_range1 = np.linspace(min(a_values_cubed)*0.98, max(a_values_cubed)*1.02, 200)
pressure_fit1 = birch_murnaghan(volume_range1, *popt1)

volume_range2 = np.linspace(min(a_values_cubed1)*0.98, max(a_values_cubed1)*1.02, 200)
pressure_fit2 = birch_murnaghan(volume_range2, *popt2)

# Plotting
plt.figure(figsize=(8, 6))

# Theoretical curve ZrC
plt.plot(pressure_ZrC_fit, volume_range_ZrC/Vo, '-', color='green', label='Birch-Murnaghan Fit ZrC')

# ML8 data and fit
plt.errorbar(pressure, a_values_cubed/popt1[0], yerr=popt1[0]*a_errors_cubed/a_values_cubed**2,
             fmt='o', color='blue', ecolor='blue', capsize=5, label='Data Zr0.4Nb0.4W0.2C')
plt.plot(pressure_fit1, volume_range1/popt1[0], '-', color='blue', label='Birch-Murnaghan fit Zr0.4Nb0.4W0.2C')

# ML3 data and fit
plt.errorbar(pressure1, a_values_cubed1/popt2[0], yerr=popt2[0]*a_errors_cubed1/a_values_cubed1**2,
             fmt='o', color="red", ecolor='red', capsize=5, label='Data ZrNbHfTaWC5')
plt.plot(pressure_fit2, volume_range2/popt2[0], '-', color='red', label='Birch-Murnaghan fit ZrNbHfTaWC5')

plt.ylabel('$V/V_0$', fontsize=12)
plt.xlabel('Pressure (GPa)', fontsize=12)
plt.xlim(0, 40)
plt.ylim(0.9,1.0)
plt.title('Birch-Murnaghan Fit: $V/V_0$ vs Pressure', fontsize=14)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
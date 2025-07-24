import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#ML8
pressure = [10.18, 18.81, 24.84, 30.29, 34.34]
a_values = np.array([4.4702, 4.4304, 4.4174, 4.3963, 4.38])
a_errors = np.array([0.0011, 0.013, 0.0027, 0.0011, 0.0054])

#ML3:
pressure1 = [3.68,9.85,15.08,20.21,25.51,31.06,35.71,]
a_values1 = np.array([4.4941,4.4652,4.4479,4.4378,4.4205,4.4071,4.3911])
a_errors1 = np.array([0.0061,0.0031,0.0035,0.0068,0.0046,0.0045,0.0075])

a_values_cubed = a_values ** 3
a_errors_cubed = 3 * (a_values ** 2) * a_errors

a_values_cubed1 = a_values1 ** 3
a_errors_cubed1 = 3 * (a_values1 ** 2) * a_errors1

#ZrC
Vo = 104.75
B = 217.7
Bo = 3.934
# Find bulk mod
def birch_murnaghan(V, V0, B0, B0_prime):
    eta = (V0 / V)
    return ((3 / 2) * B0) * (eta**(7/2) - eta**(5/2)) * (1 + (3 / 4)*(B0_prime - 4)*(eta**(2/3) - 1))
initial_guess = [a_values_cubed[0], 100, 4]
popt1, pcov1 = curve_fit(birch_murnaghan,a_values_cubed, pressure, p0=initial_guess, sigma=None, absolute_sigma=False)
V0_fitted1, B0_fitted1, B0_prime_fitted1 = popt1
errors1 = np.sqrt(np.diag(pcov1))
V0_err1, B0_err1, B0_prime_err1 = errors1


popt2, pcov2 = curve_fit(birch_murnaghan,a_values_cubed1, pressure1, p0=initial_guess, sigma=None, absolute_sigma=False)
V0_fitted2, B0_fitted2, B0_prime_fitted2 = popt2
errors2 = np.sqrt(np.diag(pcov2))
V0_err2, B0_err2, B0_prime_err2 = errors2


# Volume range for plotting
volume_range_ZrC = np.linspace(Vo * 0.1, Vo, 100)
pressure_ZrC_fit = birch_murnaghan(volume_range_ZrC, Vo, B, Bo)

print("Zr0.4Nb0.4W0.2C")
print(f"Fitted V0 = {V0_fitted1:.3f} Å³ ± {V0_err1:.3f} Å³")
print(f"Fitted Bulk modulus (B0) = {B0_fitted1:.2f} GPa ± {B0_err1:.2f} GPa")
print(f"Fitted Bulk modulus derivative (B0') = {B0_prime_fitted1:.2f} ± {B0_prime_err1:.2f}")

print("ZrNbHfTaWC5")
print(f"Fitted V0 = {V0_fitted2:.3f} Å³ ± {V0_err2:.3f} Å³")
print(f"Fitted Bulk modulus (B0) = {B0_fitted2:.2f} GPa ± {B0_err2:.2f} GPa")
print(f"Fitted Bulk modulus derivative (B0') = {B0_prime_fitted2:.2f} ± {B0_prime_err2:.2f}")
volume_range = np.linspace(min(a_values_cubed)*0.98, max(a_values_cubed)*1.02, 100)
pressure_fit = birch_murnaghan(volume_range, *popt1)

volume_range1 = np.linspace(min(a_values_cubed1)*0.98, max(a_values_cubed1)*1.02, 100)
pressure_fit1 = birch_murnaghan(volume_range1, *popt2)

plt.figure(figsize=(8, 6))
plt.plot(pressure_ZrC_fit, volume_range_ZrC, '-', color='green', label='Birch-Murnaghan Fit ZrC')

plt.errorbar(pressure,a_values_cubed, yerr=a_errors_cubed, fmt='o', color='blue', ecolor='blue', capsize=5, label='Data Zr0.4Nb0.4W0.2C')
plt.plot(pressure_fit,volume_range, '-', color='blue', label='Birch-Murnaghan fit Zr0.4Nb0.4W0.2C')
plt.errorbar(pressure1,a_values_cubed1, yerr=a_errors_cubed1, fmt='o', color ="red" ,ecolor='red', capsize=5, label='Data ZrNbHfTaWC5')
plt.plot(pressure_fit1,volume_range1, '-', color = 'red',label='Birch-Murnaghan fit ZrNbHfTaWC5')
plt.ylabel('Volume (Å³)', fontsize=12)
plt.xlabel('Pressure (GPa)', fontsize=12)
plt.xlim(0,40)
plt.ylim(60,120)
plt.title('Birch-Murnaghan Fit  Zr0.4Nb0.4W0.2C vs ZrNbHfTaWC5 Data', fontsize=14)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
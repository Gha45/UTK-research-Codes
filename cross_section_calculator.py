import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
# Define the function for the single impact model again
def single_impact(x, sigma):
    return 1 * (1 - np.exp(-sigma * x))


# Define the new data set for titanate
x_titanate = np.array([0,1.00E+11, 3.00E+11, 5.00E+11, 8.00E+11, 1.00E+12, 3.00E+12, 5.00E+12,8.00E+12])
y_titanate = np.array([0,0.044057223, 0.192689022, 0.359354668, 0.514542225, 0.5427336162, 0.987143531, 0.997650478,1])
#esds_titanate = np.array([0.00001,0.010948421, 0.031290242, 0.047833394, 0.059891973, 0.0427397984, 0.085683996, 0.100888551,0.08])
esds_titanate = np.array([0.00001,0.007213989, 0.022625563,0.044586061,0.068624131,0.04947786,0.0343273,0.005465722,0.02])

x_Ho=np.array([0,5.00E+11,8.00E+11,1.00E+12,3E12,5.00E+12,8.00E+12])
y_Ho=np.array([0,0.2561,0.2902,0.4643,0.9173,0.89202,1])
#esds_Ho=np.array([0.00001,0.0318,0.0412,0.0514,0.126,0.0868,0.0569])
esds_Ho = np.array([0.00001,0.034655443,0.038624035,0.03761539,0.020769885,0.022281311,0.02])

# look at last 5e12
# Fit the titanate data with bounds and better initial guess
popt_titanate, pcov_titanate = curve_fit(single_impact, x_titanate, y_titanate, p0=[1E-14], sigma=esds_titanate, absolute_sigma=True)
sigma_titanate = np.sqrt(np.diag(pcov_titanate))

# Fit the Ho data with bounds and better initial guess
popt_Ho, pcov_Ho = curve_fit(single_impact, x_Ho, y_Ho, p0=[4e-14], sigma=esds_Ho, absolute_sigma=True)
sigma_Ho = np.sqrt(np.diag(pcov_Ho))

print(sigma_Ho)
diameter_Ho = 2 * np.sqrt(popt_Ho[0] * 1E14 / np.pi)
diameter_err_Ho = 2 / np.sqrt(np.pi) * 0.5 * 1 / np.sqrt(popt_Ho[0] * 1E14) * sigma_Ho[0] * 1E14


dia_Ho_title = f"$SCP$"
dia_Ho = "$SCP:d_a$="+f" {round(diameter_Ho,1)} ± {round(diameter_err_Ho,1)} nm"


# Calculate diameter for titanate
diameter_titanate = 2 * np.sqrt(popt_titanate[0] * 1E14 / np.pi)
diameter_err_titanate = 2 / np.sqrt(np.pi) * 0.5 * 1 / np.sqrt(popt_titanate[0] * 1E14) * sigma_titanate[0] * 1E14


dia_ti_title = f"CCP"
dia_titanate = "$CCP:d_a$="+f" {round(diameter_titanate,1)} ± {round(diameter_err_titanate,1)} nm"

# Prepare the plot
fit_x = np.linspace(0, 5E13, 1000)
#fit_y_zirconate = single_impact(fit_x, *popt_zirconate)
fit_y_titanate = single_impact(fit_x, *popt_titanate)

fit_x_Ho = np.linspace(0, 8E12, 1500)
fit_y_Ho = single_impact(fit_x_Ho, *popt_Ho)

def plot_data_with_fits(x_data, y_data, esds_data, popt_data, color, label_prefix, shape):
    fit_x = np.linspace(0, 8E12, 1000)
    fit_y = single_impact(fit_x, *popt_data)
    plt.plot(fit_x/1E12, fit_y, label=f'{label_prefix} Fit', color=color)
    plt.errorbar(x_data/1E12, y_data, yerr=esds_data, marker=shape, ls='', color=color, label=f'{label_prefix} Data',capsize=5)

plt.figure(figsize=(10, 6))

# Plot both datasets with their respective fits
#plot_data_with_fits(x_zirconate, y_zirconate, esds_zirconate, popt_zirconate, 'blue', 'Zirconate')
plot_data_with_fits(x_titanate, y_titanate, esds_titanate, popt_titanate, 'red', 'Titanate','s')
plot_data_with_fits(x_Ho, y_Ho, esds_Ho, popt_Ho, 'black', 'Titanate','o')



# Add the text for diameters

# fontsize 16 for normal
#plt.text(0.95, 0.26, dia_ti_title, transform=plt.gca().transAxes, ha='right', va='bottom', fontsize=20, color='red')
plt.text(0.95, 0.12, dia_titanate, transform=plt.gca().transAxes, ha='right', va='bottom', fontsize=20, color='red')
#plt.text(0.95, 0.12, dia_Ho_title, transform=plt.gca().transAxes, ha='right', va='bottom', fontsize=20, color='black')
plt.text(0.95, 0.05, dia_Ho, transform=plt.gca().transAxes, ha='right', va='bottom', fontsize=20, color='black')

# Add the legend, title and labels
#plt.legend()

plt.xlabel("Fluence ($10^{12}$ ions/cm\u00b2)", fontsize=20)
plt.ylabel("Amorphous Fraction", fontsize=20)
plt.minorticks_on()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()
# Show the plot
plt.savefig("Amorphfrac.svg")
plt.show()
##

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize
import pywt
##

ACu_data = np.genfromtxt("/Users/georgeadamson/Desktop/Data/Cades data analysis/HE analysis/highentropytitanate_1e11_001.xy")

#ACu_data = np.genfromtxt("/Users/georgeadamson/Desktop/Data/Cades data analysis/HE analysis/HE_Ti2O7_5e12_001.xy")

## k(1/x^2) + b
def inversx(x, params):
    k, b = params
    return k * (1 / (x ** 2)) + b


def fit_bg_strictly_below(params, data, penalty_factor):
    residuals = data[:, 1] - inversx(data[:, 0], params)
    # Apply a very large penalty if the fit goes above the data
    penalty_above = np.sum((residuals < 0) * np.abs(residuals) * penalty_factor * 1000)
    loss = np.sum((residuals > 0) * residuals ** 2) + penalty_above
    return loss


def optimize_penalty_strict(data, initial_params, penalty_range):
    best_penalty = None
    best_loss = np.inf

    for penalty_factor in penalty_range:
        res = minimize(fit_bg_strictly_below, initial_params, args=(data, penalty_factor))
        # Check if the fit stays strictly below the data
        if res.fun < best_loss and np.all(inversx(data[:, 0], res.x) <= data[:, 1]):
            best_loss = res.fun
            best_penalty = penalty_factor

    return best_penalty


# Use this for the monstrosity definition above
penalty_range = np.linspace(100, 20000, 20)
initial_params = [0, 0]

# Optimize penalty
best_penalty_strict = optimize_penalty_strict(ACu_data, initial_params, penalty_range)

# Fit data using the best penalty factor
res_strict = minimize(fit_bg_strictly_below, initial_params, args=(ACu_data, best_penalty_strict))

# x and y values for plotting
x = ACu_data[:, 0]
y_strict = ACu_data[:, 1] - inversx(ACu_data[:, 0], res_strict.x)

# Plot results
plt.figure()
#plt.plot(ACu_data[:, 0], ACu_data[:, 1], label='data')
plt.plot(x, y_strict, lw=0.5, label='fit')
plt.plot(x, inversx(ACu_data[:, 0], res_strict.x), label='background')
plt.xlabel(r"Two-theta [A$^\circ$]", fontsize=16)
plt.ylabel(r"Intensity [counts]", fontsize=16)
plt.legend()
plt.xlim(9, 10)
plt.ylim(0, 4000)

plt.tight_layout()
plt.show()

# Print best penalty factor and final parameters
print(f"Best penalty factor: {best_penalty_strict}")
print(f"Optimized parameters: k = {res_strict.x[0]}, b = {res_strict.x[1]}")

### smoothinga data function ###
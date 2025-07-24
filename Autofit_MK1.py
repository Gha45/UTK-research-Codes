import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.signal import savgol_filter

# Load data from the uploaded file
#ACu_data = np.genfromtxt("/Users/georgeadamson/Desktop/Data/Cades data analysis/HE analysis/highentropytitanate_1e11_001.xy")
ACu_data = np.genfromtxt("/Users/georgeadamson/Desktop/Data/Cades data analysis/HE analysis/HE_Ti2O7_5e12_001.xy")

# Define the inverse x^2 model
def inversx(x, params):
    k, b = params
    return k * (1 / (x ** 2)) + b


# Define the loss function to fit below the data
def fit_bg_strictly_below(params, data, penalty_factor):
    residuals = data[:, 1] - inversx(data[:, 0], params)
    penalty_above = np.sum((residuals < 0) * np.abs(residuals) * penalty_factor * 1000)
    loss = np.sum((residuals > 0) * residuals ** 2) + penalty_above
    return loss


# Optimize penalty factor to ensure fit stays below the data
def optimize_penalty_strict(data, initial_params, penalty_range):
    best_penalty = None
    best_loss = np.inf

    for penalty_factor in penalty_range:
        res = minimize(fit_bg_strictly_below, initial_params, args=(data, penalty_factor))
        if res.fun < best_loss and np.all(inversx(data[:, 0], res.x) <= data[:, 1]):
            best_loss = res.fun
            best_penalty = penalty_factor

    return best_penalty


# Define penalty range and initial parameters
penalty_range = np.linspace(100, 20000, 20)
initial_params = [0, 0]

# Optimize penalty and fit data
best_penalty_strict = optimize_penalty_strict(ACu_data, initial_params, penalty_range)
res_strict = minimize(fit_bg_strictly_below, initial_params, args=(ACu_data, best_penalty_strict))

# Extract x and y values
x = ACu_data[:, 0]
y_strict = ACu_data[:, 1] - inversx(ACu_data[:, 0], res_strict.x)

#Savitzky-Golay filter for smoothing
window_length = 5
polyorder = 2
y_smoothed = savgol_filter(y_strict, window_length, polyorder)

# Plot original and smoothed data
plt.figure()
plt.plot(x, y_strict, lw=0.5, label='Background-subtracted data', alpha=0.5)
plt.plot(x, y_smoothed, lw=1, label='Smoothed data')
plt.xlabel(r"Two-theta [A$^\circ$]", fontsize=16)
plt.ylabel(r"Intensity [counts]", fontsize=16)
plt.legend()
plt.xlim(4.5, 11)
plt.ylim(1, 200)
plt.tight_layout()
plt.show()

# Print the best penalty factor and optimized parameters
best_penalty_strict, res_strict.x


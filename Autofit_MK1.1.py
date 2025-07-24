import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize
from scipy.signal import find_peaks_cwt

# ---------------------------------------------------------
# Step 1: Load Data and Subtract Background
# ---------------------------------------------------------
# Replace with the actual file path for your data.
ACu_data = np.genfromtxt("/Users/georgeadamson/Desktop/Data/Cades data analysis/HE analysis/HE_Ti2O7_5e12_001.xy")

# Background model: k * (1/x^2) + b
def inversx(x, params):
    k, b = params
    return k * (1 / (x ** 2)) + b

# Loss function that penalizes fits which exceed the data
def fit_bg_strictly_below(params, data, penalty_factor):
    residuals = data[:, 1] - inversx(data[:, 0], params)
    penalty_above = np.sum((residuals < 0) * np.abs(residuals) * penalty_factor * 1000)
    loss = np.sum((residuals > 0) * residuals ** 2) + penalty_above
    return loss

# Optimize the penalty factor over a range of values
def optimize_penalty_strict(data, initial_params, penalty_range):
    best_penalty = None
    best_loss = np.inf
    for penalty_factor in penalty_range:
        res = minimize(fit_bg_strictly_below, initial_params, args=(data, penalty_factor))
        if res.fun < best_loss and np.all(inversx(data[:, 0], res.x) <= data[:, 1]):
            best_loss = res.fun
            best_penalty = penalty_factor
    return best_penalty

penalty_range = np.linspace(100, 20000, 20)
initial_params = [0, 0]
best_penalty_strict = optimize_penalty_strict(ACu_data, initial_params, penalty_range)
res_strict = minimize(fit_bg_strictly_below, initial_params, args=(ACu_data, best_penalty_strict))

# Get x (two-theta) and subtract the fitted background
x = ACu_data[:, 0]
y_bg = inversx(x, res_strict.x)
y_strict = ACu_data[:, 1] - y_bg

# Plot raw data, background, and background‐subtracted data.
plt.figure(figsize=(10, 4))
plt.plot(x, ACu_data[:, 1], label="Raw Data", alpha=0.6)
plt.plot(x, y_bg, label="Fitted Background", linestyle="--")
plt.plot(x, y_strict, label="Background Subtracted", lw=0.5)
plt.xlabel(r"Two-theta [A$^\circ$]", fontsize=14)
plt.ylabel("Intensity [counts]", fontsize=14)
plt.legend()
plt.xlim(9, 10)
plt.ylim(0, 4000)
plt.tight_layout()
plt.show()

print("Best penalty factor:", best_penalty_strict)
print("Optimized background parameters: k = {:.4f}, b = {:.4f}".format(res_strict.x[0], res_strict.x[1]))

# ---------------------------------------------------------
# Step 2: Wavelet-based Peak Detection
# ---------------------------------------------------------
# Use a range of widths (scales) for the wavelet transform.
widths = np.arange(1, 20)
detected_peak_indices = find_peaks_cwt(y_strict, widths)

# Get detected peak positions and amplitudes.
detected_peak_positions = x[detected_peak_indices]
detected_peak_amplitudes = y_strict[detected_peak_indices]

print("Number of peaks detected (via wavelets):", len(detected_peak_positions))
print("Detected peak positions:", detected_peak_positions)

# ---------------------------------------------------------
# Step 3: Prompt for the Number of Peaks for Each Model
# ---------------------------------------------------------
# Total combined model: some peaks modeled as Gaussian and others as Pseudovoigt.
num_gaussian = int(input("Enter the number of Gaussian peaks you expect: "))
num_pseudovoigt = int(input("Enter the number of Pseudovoigt peaks you expect: "))
total_peaks_expected = num_gaussian + num_pseudovoigt

# From the detected peaks, select the top total_peaks_expected peaks by amplitude.
if len(detected_peak_positions) >= total_peaks_expected:
    sort_order = np.argsort(detected_peak_amplitudes)[::-1]  # descending order by amplitude
    selected_indices = [detected_peak_indices[i] for i in sort_order[:total_peaks_expected]]
else:
    selected_indices = detected_peak_indices

# Sort the selected indices in increasing order of x.
selected_indices = np.sort(selected_indices)
selected_positions = x[selected_indices]
selected_amplitudes = np.array([y_strict[(np.abs(x - pos)).argmin()] for pos in selected_positions])

print("Total selected peak positions:", selected_positions)

# Assign the first num_gaussian peaks to the Gaussian model,
# and the remaining to the Pseudovoigt model.
selected_positions_gauss = selected_positions[:num_gaussian]
selected_amplitudes_gauss = selected_amplitudes[:num_gaussian]

selected_positions_pseudo = selected_positions[num_gaussian:]
selected_amplitudes_pseudo = selected_amplitudes[num_gaussian:]

print("Gaussian peak positions:", selected_positions_gauss)
print("Pseudovoigt peak positions:", selected_positions_pseudo)

# ---------------------------------------------------------
# Step 4: Define Peak Functions
# ---------------------------------------------------------
def gaussian(x, A, mu, sigma):
    """Gaussian peak function."""
    return A * np.exp(-((x - mu)**2) / (2 * sigma**2))

def lorentzian(x, A, mu, gamma):
    """Lorentzian peak function."""
    return A * (gamma**2) / ((x - mu)**2 + gamma**2)

def pseudovoigt(x, A, mu, sigma, eta):
    """
    Pseudovoigt peak function: a linear combination of Gaussian and Lorentzian.
    'eta' (0 <= eta <= 1) weights the Lorentzian component.
    """
    return eta * lorentzian(x, A, mu, sigma) + (1 - eta) * gaussian(x, A, mu, sigma)

# ---------------------------------------------------------
# Step 5: Define Composite Functions for Each Model Component
# ---------------------------------------------------------
def multi_gaussian(x, *params):
    """
    Sum of multiple Gaussian peaks.
    For each peak, parameters: [A, mu, sigma]
    """
    n_peaks = len(params) // 3
    y = np.zeros_like(x)
    for i in range(n_peaks):
        A = params[3*i]
        mu = params[3*i + 1]
        sigma = params[3*i + 2]
        y += gaussian(x, A, mu, sigma)
    return y

def multi_pseudovoigt(x, *params):
    """
    Sum of multiple Pseudovoigt peaks.
    For each peak, parameters: [A, mu, sigma, eta]
    """
    n_peaks = len(params) // 4
    y = np.zeros_like(x)
    for i in range(n_peaks):
        A = params[4*i]
        mu = params[4*i + 1]
        sigma = params[4*i + 2]
        eta = params[4*i + 3]
        y += pseudovoigt(x, A, mu, sigma, eta)
    return y

# ---------------------------------------------------------
# Step 6: Define Combined Model Function
# ---------------------------------------------------------
def combined_model(x, *params):
    """
    Combined model: sum of Gaussian and Pseudovoigt peaks.
    The parameters are arranged as:
      - First 3*num_gaussian parameters for Gaussian peaks.
      - Next 4*num_pseudovoigt parameters for Pseudovoigt peaks.
    """
    gauss_params = params[:3*num_gaussian]
    pseudo_params = params[3*num_gaussian:]
    return multi_gaussian(x, *gauss_params) + multi_pseudovoigt(x, *pseudo_params)

# ---------------------------------------------------------
# Step 7: Build Initial Guesses and Bounds for the Combined Fit
# ---------------------------------------------------------
# For Gaussian peaks: [Amplitude, mu, sigma]
initial_guess_gauss = []
lower_bounds_gauss = []
upper_bounds_gauss = []
for pos, amp in zip(selected_positions_gauss, selected_amplitudes_gauss):
    initial_guess_gauss.extend([amp, pos, 0.05])
    lower_bounds_gauss.extend([0, pos - 0.1, 0.001])
    upper_bounds_gauss.extend([np.inf, pos + 0.1, 1.0])

# For Pseudovoigt peaks: [Amplitude, mu, sigma, eta]
initial_guess_pseudo = []
lower_bounds_pseudo = []
upper_bounds_pseudo = []
for pos, amp in zip(selected_positions_pseudo, selected_amplitudes_pseudo):
    initial_guess_pseudo.extend([amp, pos, 0.05, 0.5])
    lower_bounds_pseudo.extend([0, pos - 0.1, 0.001, 0])
    upper_bounds_pseudo.extend([np.inf, pos + 0.1, 1.0, 1])

# Combine the initial guesses and bounds
initial_guess_total = initial_guess_gauss + initial_guess_pseudo
lower_bounds_total = lower_bounds_gauss + lower_bounds_pseudo
upper_bounds_total = upper_bounds_gauss + upper_bounds_pseudo

# ---------------------------------------------------------
# Step 8: Fit the Combined Model Using curve_fit
# ---------------------------------------------------------
maxfev_val = 100000
popt_total, pcov_total = curve_fit(
    combined_model, x, y_strict,
    p0=initial_guess_total,
    bounds=(lower_bounds_total, upper_bounds_total),
    maxfev=maxfev_val
)

# ---------------------------------------------------------
# Step 9: Plot the Results and Print Fitted Parameters
# ---------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.plot(x, y_strict, label="Background Subtracted Data", color="black")
plt.plot(x, combined_model(x, *popt_total), label="Combined Fit", linestyle="--")
plt.xlabel(r"Two-theta [A$^\circ$]", fontsize=14)
plt.ylabel("Intensity [counts]", fontsize=14)
plt.legend()
plt.tight_layout()
plt.show()

# Print fitted parameters for Gaussian peaks.
print("Fitted Gaussian Peaks:")
for i in range(num_gaussian):
    A, mu, sigma = popt_total[3*i:3*i+3]
    print("Gaussian Peak {}: Amplitude = {:.2f}, Position = {:.4f}, Sigma = {:.4f}".format(i+1, A, mu, sigma))

# Print fitted parameters for Pseudovoigt peaks.
print("\nFitted Pseudovoigt Peaks:")
for i in range(num_pseudovoigt):
    idx = 3*num_gaussian + 4*i
    A, mu, sigma, eta = popt_total[idx:idx+4]
    print("Pseudovoigt Peak {}: Amplitude = {:.2f}, Position = {:.4f}, Sigma = {:.4f}, eta = {:.4f}".format(i+1, A, mu, sigma, eta))

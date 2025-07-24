import os
import numpy as np
import matplotlib.pyplot as plt

# Specify the folder path containing the .txt files
folder = '/Users/georgeadamson/Desktop/UTK python programs for reaserch/Data/ErHo2TI2O7redraman'
z = 0
# cost function
def select_edge_points(x_data, y_data, left_start=0, left_end=50, left_points=3, right_start=1150, right_points=5):
    """
    Select points from the leftmost side between specified indices and a specified number of points
    evenly spaced from a starting index on the right side to the end of the dataset.
    """
    if left_end >= len(y_data) or left_end < left_start:
        raise ValueError("Invalid range for left point selection.")
    if right_start >= len(y_data) or len(y_data) - right_start < right_points:
        raise ValueError("Right start index is too high or not enough points to select from.")

    # Calculate even spacing for left-hand points
    left_indices = np.linspace(left_start, left_end, left_points, dtype=int)

    # Calculate even spacing for right-hand points
    total_right_points = len(y_data) - right_start
    step = total_right_points // right_points
    right_indices = np.arange(right_start, len(y_data), step)[:right_points]

    # Combine indices from both sides
    indices = np.unique(np.concatenate((left_indices, right_indices)))
    return indices, x_data[indices], y_data[indices]

def fit_constrained_polynomial(x_data, y_data, indices, degree=2, max_iterations=20):
    """
    Fit a polynomial to the provided points, ensuring it does not exceed the actual data points.
    """
    x_points = x_data[indices]
    y_points = y_data[indices]
    for iteration in range(max_iterations):
        coefficients = np.polyfit(x_points, y_points, degree)
        polynomial = np.poly1d(coefficients)
        polynomial_fit = polynomial(x_data)
        if np.all(polynomial_fit <= y_data):
            break
        # Slightly reduce y-values to ensure the polynomial stays below the data points
        y_points -= 0.01 * np.min(y_points)
    return polynomial

# Process files in the specified folder
for filename in os.listdir(folder):
    z += 1
    if filename.endswith(".txt"):
        file_path = os.path.join(folder, filename)
        data_array = np.genfromtxt(file_path, delimiter='\t', encoding='latin-1')
        if data_array.ndim > 1 and data_array.shape[1] > 1:
            x_data = np.arange(len(data_array[:, 1]))
            y_data = data_array[:, 1]
            indices, x_points, y_points = select_edge_points(x_data, y_data, left_start=0, left_end=50, left_points=2, right_start=1150, right_points=10)
            polynomial = fit_constrained_polynomial(x_data, y_data, indices)
            polynomial_fit = polynomial(x_data)
            # Plotting
            corrected_data = y_data - polynomial_fit

            # Plotting results for verification
            plt.figure(figsize=(12, 6))
            plt.subplot(2, 1, 1)
            plt.plot(x_data, y_data, label='Original Data')
            plt.plot(x_data, polynomial_fit, label='Polynomial Background', linestyle='--')
            plt.scatter(x_points, y_points, color='red', label='Selected Points', s=50, zorder=5)
            plt.legend()
            plt.title(f'Original Data and Polynomial Background for {filename}')

            plt.subplot(2, 1, 2)
            plt.plot(x_data, corrected_data, label='Corrected Data', color='green')
            plt.legend()
            plt.title(f'Corrected Data After Background Subtraction for {filename}')

            plt.xlabel('Data Points')
            plt.ylabel('Intensity')
            plt.tight_layout()
            plt.show()

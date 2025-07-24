import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
from scipy.optimize import curve_fit
from scipy.integrate import simpson
from numpy import trapz
# Specify the folder path containing the .txt files
#folder =  "/Users/georgeadamson/Desktop/Data/Raman_data/2024 raman data/Removal_of_copies_He_pyrochlore_2024 copy/YbErHoDyTb2Ti2O7"
folder1 = "/Users/georgeadamson/Desktop/Data/Raman_data/Raman_for_CCO_paper/Ho"
folder2 ='/Users/georgeadamson/Desktop/Data/Raman_data/Raman_for_CCO_paper/HE'
#folder ='/Users/georgeadamson/Desktop/Data/Raman_data/2024 RAMAN /DyGdEuSmNd2Ti2O7'
# Specify the column ranges for calculating max and min values from the second column
# fitting coreseponds to the differnce to the maximal and minimal slope of the exponential curve, using the errors
max_column_range = [150, 350]
min_column_range = [315, 500]

# Initialize dictionaries to store data
pattern_averages = {}
pattern_max_values = {}
pattern_min_values = {}
pattern_data = {}

for filename in os.listdir(folder):
    if filename.endswith(".txt"):
        file_path = os.path.join(folder, filename)
        data_array = np.genfromtxt(file_path, delimiter='\t', encoding='latin-1')

        if data_array.ndim > 1 and data_array.shape[1] > 1:
            wavenumber_data = data_array[:, 0]
            column_data = data_array[:, 1]

            # Normalize the data by setting the area under the curve to 1 for each
            area = trapz(column_data, dx = len(column_data))
            #print("area =", area)
            normalized_column = column_data/area
            #normalized_column = column_data
            #normed_area  = trapz(normalized_column, dx = len(column_data))
            #print('norm area:',normed_area)

            # Calculate max and min values within specified ranges
            max_value = np.max(normalized_column[max_column_range[0]:max_column_range[1]])
            min_value = np.min(normalized_column[min_column_range[0]:min_column_range[1]])

            # Determine pattern from filename
            pattern = 0.0 if 'pris' in filename else float(f"{m.group(1)}E{m.group(2)}") if (
                m := re.search(r'(\d+)[Ee](\d+)', filename)) else None
            if pattern is not None:
                pattern_averages.setdefault(pattern, []).append(max_value / min_value if min_value != 0 else None)
                pattern_max_values.setdefault(pattern, []).append(max_value)
                pattern_min_values.setdefault(pattern, []).append(min_value)
                pattern_data.setdefault(pattern, []).append((wavenumber_data, normalized_column))

# Calculate statistics for each pattern
data = []
for pattern, ratios in pattern_averages.items():
    max_values = pattern_max_values[pattern]
    min_values = pattern_min_values[pattern]

    max_mean = np.mean(max_values)
    min_mean = np.mean(min_values)
    max_std = np.std(max_values)
    min_std = np.std(min_values)

    # Custom standard deviation calculation
    std_dev = np.sqrt(((max_std / max_mean) ** 2) + (((-max_mean * min_std) / (min_mean ** 2)) ** 2))

    # Check if std_dev is zero and adjust if necessary
    if std_dev == 0:
        average_ratio = np.mean(ratios)
        std_dev = 0.05 * average_ratio


    #print(std_dev)

    # Storing results in a dictionary
    data.append({
        'Pattern': pattern,
        'Average Ratio': np.mean(ratios),
        'Standard Deviation': std_dev  # Uncomment normalization if required
    })
# Create DataFrame, sort by 'Pattern', and reset the index
df = pd.DataFrame(data).sort_values('Pattern').reset_index(drop=True)

# Extract lists from DataFrame
Pattern_list = df['Pattern'].tolist()
AverageRatio_list = df['Average Ratio'].tolist()
StandardDeviation_list = df['Standard Deviation'].tolist()


# Define the model functions
def gen_single_impact(start):
    def single_imp(x, sigma, b):
        return start + b * (1 - np.exp(-sigma * x))

    return single_imp


def calculate_diameter_Ho(gen_single_impact, x_Ho, y_Ho, esds_Ho):
    initial_guess = [3E-12, -3]
    # make line with the minimal slope, line with maximal slope and take the distance between the two
    # print(x_data)
    # print(y_data)
    y_min = [y - e for y, e in zip(y_Ho, esds_Ho)]
    # print(y_min)
    y_max = [y + e for y, e in zip(y_Ho, esds_Ho)]
    # print(y_max)
    # for line 1
    y_Ho_updated_1 = y_Ho.copy()
    y_Ho_updated_1[0] = y_min[0]
    y_Ho_updated_1[1] = y_max[1]
    y_Ho_updated_1[2] = y_max[2]

    # for line 2
    y_Ho_updated_2 = y_Ho.copy()
    y_Ho_updated_2[0] = y_max[0]
    y_Ho_updated_2[1] = y_min[1]
    y_Ho_updated_2[2] = y_min[2]

    # for the actual value
    popt_Ho_1, pcov_Ho_1 = curve_fit(gen_single_impact(y_Ho[0]), x_Ho, y_Ho, p0=initial_guess, sigma=esds_Ho)
    popt_Ho_2, pcov_Ho_2 = curve_fit(gen_single_impact(y_Ho_updated_1[0]), x_Ho, y_Ho_updated_1, p0=initial_guess)
    popt_Ho_3, pcov_Ho_3 = curve_fit(gen_single_impact(y_Ho_updated_2[0]), x_Ho, y_Ho_updated_2, p0=initial_guess)
    sigma_Ho_1 = np.sqrt(np.diag(pcov_Ho_1))
    sigma_Ho_2 = np.sqrt(np.diag(pcov_Ho_2))
    sigma_Ho_3 = np.sqrt(np.diag(pcov_Ho_2))
    diameter_Ho_1 = 2 * np.sqrt(popt_Ho_1[0] * 1E14 / np.pi)
    diameter_Ho_2 = 2 * np.sqrt(popt_Ho_2[0] * 1E14 / np.pi)
    diameter_Ho_3 = 2 * np.sqrt(popt_Ho_3[0] * 1E14 / np.pi)

    diameter_err_Ho = (abs(diameter_Ho_3 - diameter_Ho_2))/2
    #diameter_err_Ho = 1 / np.sqrt(np.pi) * 1 / np.sqrt(popt_Ho[0] * 1E14) * (sigma_Ho[0] * 1E14)

    #print("Optimized parameters:", popt_Ho)
    #print("Parameter errors:", sigma_Ho)
    #print("Diameter:", diameter_Ho)
    #print('Error:', diameter_err_Ho)
    dia_Ho = "$D$=" + f" {round(diameter_Ho_1, 1)} ± {round(diameter_err_Ho, 1)} nm"
    #print(dia_Ho)
    return popt_Ho_1,popt_Ho_2,popt_Ho_3 , y_Ho_updated_1 ,  y_Ho_updated_2, sigma_Ho_1, diameter_Ho_1, diameter_err_Ho, dia_Ho


def plot_data_with_fits(x_data, y_data, esds_data, popt_1,popt_2,popt_3,y_Ho_updated_1, y_Ho_updated_2 , color, label_prefix, dia_Ho):
    x_data = np.array(x_data)
    y_data = np.array(y_data)
    esds_data = np.array(esds_data)


    fit_x = np.linspace(0, 1E13, 1000)
    fit_y = gen_single_impact(y_data[0])(fit_x, *popt_1)
    fit_y1 = gen_single_impact(y_Ho_updated_1[0])(fit_x, *popt_2)
    fit_y2 = gen_single_impact(y_Ho_updated_2[0])(fit_x, *popt_3)


    plt.plot(fit_x / 1E12, fit_y, label=f'{label_prefix} Fit', color=color)
    #plt.plot(fit_x / 1E12, fit_y1, label=f'{label_prefix} Fit', color='green')
    #plt.plot(fit_x / 1E12, fit_y2, label=f'{label_prefix} Fit', color='purple')
    plt.errorbar(x_data / 1E12, y_data, yerr=esds_data, marker='s', ls='', color="blue", label=f'{label_prefix} Data')
    plt.xlim(-0.1, 10)
    plt.xlabel("Fluence ($10^{12}$ ions/cm\u00b2)", fontsize=14)
    plt.ylabel("Intenisty Ratio", fontsize=14)

    #"Ho$_2$Ti$_2$O$_7$"
    # "(YbTbGdErDy)$_2$Ti$_2$O$_7$"
    plt.text(0.95, 0.15, "Ho$_2$Ti$_2$O$_7$", transform=plt.gca().transAxes, ha='right', va='bottom', fontsize=14,
             color='red')
    plt.text(0.95, 0.1, dia_Ho, transform=plt.gca().transAxes, ha='right', va='bottom', fontsize=14, color='red')
    plt.savefig('DyGdEuSmNd2Ti2O7_GSI2024_Dia_Raman.svg')
    plt.show()


def plot_stacked_patterns_per_fluence(pattern_data):
    print(pattern_data)
    for pattern, data_list in pattern_data.items():
        plt.figure(figsize=(10, 10))
        for wavenumber_data, normalized_data in data_list:
            plt.plot(wavenumber_data, normalized_data, label=f'Pattern {format_e(pattern)}')
        plt.xlabel('Wavenumber')
        plt.ylabel('Normalized Intensity Ratio')
        plt.title(f'Stacked Plot for Pattern {format_e(pattern)}')
        plt.legend(loc='upper right')

        plt.show()
def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

# Perform calculations and plot
popt_Ho = calculate_diameter_Ho(gen_single_impact, Pattern_list, AverageRatio_list, StandardDeviation_list)
plot_data_with_fits(Pattern_list, AverageRatio_list, StandardDeviation_list, popt_Ho[0],popt_Ho[1],popt_Ho[2],popt_Ho[3],popt_Ho[4], 'red', 'test', popt_Ho[8])

#plot_stacked_patterns_per_fluence(pattern_data)

# Create a DataFrame with the necessary values
output_data = {
    "Fluence (Pattern)": Pattern_list,
    "Intensity Ratio (Y)": AverageRatio_list,
    "Error (Standard Deviation)": StandardDeviation_list,
    "Diameter (nm)": [popt_Ho[6]] * len(Pattern_list),  # Diameter is constant for all
    "Diameter Error (nm)": [popt_Ho[7]] * len(Pattern_list)  # Same for error
}

df_output = pd.DataFrame(output_data)

# Define CSV file name
csv_filename = "Raman_analysis_output.csv"

# Save DataFrame to CSV
df_output.to_csv(csv_filename, index=False)

print(f"CSV file saved: {csv_filename}")


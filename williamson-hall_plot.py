import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import re
import math
from scipy.odr import Model, Data, ODR

def linear(params, x):
    return params[0] + params[1] * x
# General flow:
# 1. Make Model from your selected function
# 2. Make RealData using your uncertainties
# 3. Run ODR on the Model and provide some starting values
# 4. Print fit output and use refined parameters to plot your fit

def linear_model(x, a, b):
    return a * x + b
def quadratic_model(x, a, b, c):
    return a * x**2 + b * x + c

def read_complex_data(file_path):
    # Open and read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Define a pattern to match "NUMBER +/- NUMBER"
    pattern = re.compile(r'([\d\.]+) \+/- ([\d\.e\+-]+)')

    # Parse the data
    data = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) > 5 and "Gaussian" not in parts:  # Filter out rows with "Gaussian"
            # Extract the fixed number of fields
            peak_type = parts[1]  # This assumes the first part is a marker like '%_1'
            center = parts[2]
            height = parts[3]
            area = parts[4]
            fwhm = parts[5]
            # Handle the parameters part
            parameters = ' '.join(parts[6:])
            matches = pattern.findall(parameters)
            # Append parsed data and ensure we have 4 matches, else fill with None
            data.append([peak_type, center, height, area, fwhm] + [m for match in matches for m in match] + [None] * (
                        8 - 2 * len(matches)))

    # Create DataFrame with extended columns
    columns = ['PeakType', 'Center', 'Height', 'Area', 'FWHM'] + [f'Param_{i}{suffix}' for i in range(1, 5) for suffix
                                                                  in ('_Value', '_Error')]
    df = pd.DataFrame(data, columns=columns)
    return df


file_path ='//Users/georgeadamson/Desktop/UTK python programs for reaserch/Williamson_Hall_peaks/Pristine GSI 2020/Ho_Ti/Ho2Ti2O7_pristine_20210617-121226_20e625_0001_mean_tth.peaks'
#file_path = '/Users/georgeadamson/Desktop/UTK python programs for reaserch/Williamson_Hall_peaks/GSI2024-NSLS-2/Peaks/Titanate/LuYbEr2Ti2O7/xrd_LuYbEr2Ti2O71400CSafin_20240318-164214_25c7e2_primary-1_mean_tth.peaks'
# Use the function to read the data
df = read_complex_data(file_path)

df = df.dropna()
threshold_value = 0.35

df['Center'] =df.loc[1:, 'Center']
df['Center'] = np.radians(df['Center'].astype(float))/2
df['Center_error'] =np.radians(df['Param_2_Error'].astype(float))
df['Param_3_Value'] = np.radians(df['Param_3_Value'].astype(float)*2)
df['FWHM_error'] = np.radians((df['Param_3_Error'].astype(float)) * 2)


df['dBeta'] = np.sqrt((df['Center_error'])**2+(df['FWHM_error'])**2)

df['x-axis'] = 4 * np.sin(df['Center'])
df = df[df['x-axis'] <= threshold_value]

#GSI 2020
#intermediate_calculation = ((df['Param_3_Value'].astype(float))**2) - (((0.0007099294692323462 * (df['Center'].astype(float))**2) + ((-1.3177335197514796e-05) * df['Center'].astype(float)) + 0.0003974675988717483)**2)
#df['y-axis'] = np.sqrt(intermediate_calculation) * np.cos(df['Center'].astype(float))
#GSI 2024
#intermediate_calculation = ((df['FWHM'].astype(float))**2) - (((0.00135120 * (df['Center'].astype(float))**2) - ((0.00004910) * df['Center'].astype(float)) + 0.0004874)**2)
#df['y-axis'] = np.sqrt(intermediate_calculation) * np.cos(df['Center'].astype(float))
#df['y-axis'] = df['FWHM'].astype(float) * np.cos(df['Center'])
#2021

intermediate_calculation = ((df['Param_3_Value'].astype(float))**2) - (((0.000811175888 * (df['Center'].astype(float))**2) + (( 0.000022872696) * df['Center'].astype(float)) + 0.000308591696)**2)
df['y-axis'] = np.sqrt(intermediate_calculation) * np.cos(df['Center'].astype(float))
def remove_outliers(df, threshold=0.45):
    """
    Removes rows in the DataFrame where the percentage difference between
    consecutive y-axis values exceeds a given threshold.

    Parameters:
    df (pd.DataFrame): The DataFrame containing 'y-axis' values.
    threshold (float): The threshold percentage for detecting outliers (default is 10%).

    Returns:
    pd.DataFrame: The DataFrame with outliers removed.
    """
    # Initialize a list to store indices of rows to keep
    indices_to_keep = [0]  # Always keep the first row

    # Loop through the DataFrame starting from the second row
    for i in range(1, len(df)):
        # Calculate the percentage difference from the previous row
        prev_y = df['y-axis'].iloc[i - 1]
        current_y = df['y-axis'].iloc[i]
        percentage_diff = abs((current_y - prev_y) / prev_y)

        # If the percentage difference is less than the threshold, keep the row
        if percentage_diff <= threshold:
            indices_to_keep.append(i)

    # Return the DataFrame with only the rows that are not considered outliers
    return df.iloc[indices_to_keep].reset_index(drop=True)
df = df.sort_values(by='x-axis').reset_index(drop=True)
df = remove_outliers(df)


#integral error: -beta* sintheta * dtheta * costheta dbeta

#df['y-error'] = np.radians(df['Param_3_Error'].astype(float))   # Assuming this is already in appropriate units
df['y-error'] = np.cos(df['Center']) * df['dBeta']
x = df['x-axis'][1:].reset_index(drop=True)
y = df['y-axis'][1:].reset_index(drop=True)
sy = df['y-error'][1:].reset_index(drop=True)

#Max slope
y_min = [y - e for y, e in zip(y, sy)]
y_max = [y + e for y, e in zip(y, sy)]

y_updated_1 = y.copy()
y_updated_1.iloc[0] = y_max[0]
#y_updated_1.iloc[1] = y_max[1]
#y_updated_1.iloc[2] = y_max[2]

y_updated_1.iloc[-1] = y_min[-1]
#y_updated_1.iloc[-2] = y_min[-2]
#y_updated_1.iloc[-3] = y_min[-3]

# for line 2


# linear
popt, pcov = curve_fit(linear_model, x, y,sigma=sy)
intercept_error_fit = np.sqrt(pcov[1, 1])


popt1, pcov1 = curve_fit(linear_model, x, y_max)
popt2, pcov2 = curve_fit(linear_model, x, y_min)

y_fit_curve = linear_model(x, *popt)
y_fit_curve1 = linear_model(x, *popt1)
y_fit_curve2 = linear_model(x, *popt2)
slope, intercept = popt
per = (intercept_error_fit/intercept)/10


print(f"Slope: {slope}, Intercept: {intercept}")
equation = f"y = {slope}x + {intercept}"
print("Equation of the fitted line:", equation)

err1 = (0.9*0.1665)/popt1[1]
err2 = (0.9*0.1665)/popt2[1]
intercept_err = (abs(err2 - err1)/2)

# quadratic
#popt, pcov = curve_fit(quadratic_model, x, y)
#y_fit_curve = quadratic_model(x, *popt)
#a,b,c = popt
#print(a,b,c)

#equation = f"$y = {a:.8f}x^2 + {b:.8f}x + {c:.8f}$"
#print("Equation of the fitted line:", equation)

grainsize = (0.9*0.1665)/intercept
print("err",(grainsize * per))
grain_text = f"grainsize(nm):{grainsize/10:.1f}" #±{intercept_err/10:.1f}
strain = slope/((10**-6))
strain_text = f"microstrain:{strain:.1f}"
print(strain_text)
#grain_text = f"grainsize(nm):{grainsize/10:.1f}"

# Optionally, you can plot the result
import matplotlib.pyplot as plt


plt.figure(figsize=(10, 6))
plt.tight_layout()
plt.plot(x, y_fit_curve)
#plt.plot(x, y_fit_curve1)
#plt.plot(x, y_fit_curve2)

plt.text(0.05, 0.95, grain_text, transform=plt.gca().transAxes,
         verticalalignment='top', fontsize=20, bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white'))
plt.errorbar(x, y, yerr=sy, linestyle='None', marker='o', label='Data with Errors')
plt.text(0.05, 0.85, strain_text, transform=plt.gca().transAxes,
         verticalalignment='top', fontsize=20, bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white'))

plt.xlabel('(4*sin(Theta))',fontsize = 20)
plt.xticks(fontsize = 20)
plt.ylabel('B*cos(Theta)',fontsize = 20)
plt.yticks(fontsize = 20)
plt.title("Ho2Ti2O7")

plt.show()
#df.to_csv('testing.csv', sep=',')
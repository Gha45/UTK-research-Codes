import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import re
import math
from scipy.odr import Model, Data, ODR
import matplotlib.pyplot as plt


def linear_model(x, a, b):
    return a * x + b

def quadratic_model(x, a, b, c):
    return a * x**2 + b * x + c

def read_complex_data(file_path):
    """
    Reads a .peaks file, extracts columns [PeakType, Center, Height, Area, FWHM, ...].
    Ensures that we have columns like Param_1_Value, Param_1_Error, etc.
    Returns a DataFrame.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    pattern = re.compile(r'([\d\.]+) \+/- ([\d\.e\+-]+)')

    data = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) > 5 and "Gaussian" not in parts:  # Filter out rows containing "Gaussian"
            peak_type = parts[1]  # E.g. '%_1'
            center    = parts[2]
            height    = parts[3]
            area      = parts[4]
            fwhm      = parts[5]
            parameters = ' '.join(parts[6:])
            matches    = pattern.findall(parameters)
            # Each match is like [('0.123','0.001'),('4.56','0.89'),...]
            # Flatten into [v1,e1,v2,e2, ...]
            flattened = []
            for match in matches:
                flattened.extend(match)
            # Ensure we have exactly 8 columns for param values, fill with None if missing
            while len(flattened) < 8:
                flattened.append(None)
            row = [peak_type, center, height, area, fwhm] + flattened[:8]
            data.append(row)

    columns = ['PeakType', 'Center', 'Height', 'Area', 'FWHM'] + [
        f'Param_{i}{suffix}' for i in range(1, 5) for suffix in ('_Value', '_Error')
    ]
    df = pd.DataFrame(data, columns=columns)
    return df

def remove_outliers(df, threshold=1):
    """
    Removes rows in the DataFrame where the percentage difference between
    consecutive y-axis values exceeds `threshold`.
    """
    indices_to_keep = [0]  # Always keep the first row
    for i in range(1, len(df)):
        prev_y = df['y-axis'].iloc[i - 1]
        curr_y = df['y-axis'].iloc[i]
        if prev_y == 0:
            continue
        perc_diff = abs((curr_y - prev_y) / prev_y)
        if perc_diff <= threshold:
            indices_to_keep.append(i)
    return df.iloc[indices_to_keep].reset_index(drop=True)

def process_peaks_file(file_path, label, y_calc_fn):
    """
    Reads file_path, does all the pre-processing, uses `y_calc_fn(df)`
    to define the final 'y-axis' column, then does a linear fit.
    Returns a dictionary with {x, y, sy, y_fit, grain_text, strain_text, label}.
    """

    # ---------- Read & drop blanks -----------
    df = read_complex_data(file_path).dropna()

    # ---------- Convert angles -----------
    # We skip the first row for 'Center' in your script (df.loc[1:, 'Center']), but let's do it consistently
    df['Center']       = np.radians(df.loc[1:, 'Center'].astype(float))/2
    df['Center_error'] = np.radians(df['Param_2_Error'].astype(float))
    # Param_3_Value is typically the FWHM in degrees, so we convert it to radians
    df['Param_3_Value'] = np.radians(df['Param_3_Value'].astype(float) * 2)
    # If you also need the 'FWHM' column specifically in the broadening formula:
    # make sure it's a float in radians as well
    df['FWHM_error']    = np.radians(df['Param_3_Error'].astype(float) * 2)

    # ---------- Calculate dBeta -----------
    df['dBeta'] = np.sqrt(df['Center_error']**2 + df['FWHM_error']**2)

    # ---------- X-axis: 4*sin(theta) -----------
    df['x-axis'] = 4 * np.sin(df['Center'])

    # ---------- Threshold -----------
    threshold_value = 0.30
    df = df[df['x-axis'] <= threshold_value]

    # ---------- y_calc_fn for instrument broadening -----------
    df['y-axis'] = y_calc_fn(df)

    # ---------- Remove outliers & reindex -----------
    df = df.sort_values('x-axis').reset_index(drop=True)
    df = remove_outliers(df)

    # ---------- Propagate error in y (approx. = cos(theta)*dBeta) -----------
    df['y-error'] = np.cos(df['Center']) * df['dBeta']

    # ---------- Trim out the first row after reindex if needed -----------
    # The original script was skipping df[0], so let's do that here too
    x  = df['x-axis'][1:].reset_index(drop=True).astype(float)
    y  = df['y-axis'][1:].reset_index(drop=True).astype(float)
    sy = df['y-error'][1:].reset_index(drop=True).astype(float)

    # ---------- Fit the linear model: y = slope*x + intercept -----------
    popt, pcov = curve_fit(linear_model, x, y, sigma=sy, absolute_sigma=True)
    slope, intercept = popt
    y_fit = linear_model(x, *popt)
    slope_e, intercept_e = np.polyfit(x, y, 1)

    # ---------- Some reporting -----------
    # (grain size in nm) = (K * lambda) / (intercept * cos(theta)?)
    # Based on your usage: K=0.9, lambda=0.1665, dividing by intercept*10 => nm
    grain_size_nm = (0.9*0.1665)/(intercept*10)
    microstrain = slope / (1e-6)  # slope * 1e6 => microstrain

    grain_text  = f"grainsize: {grain_size_nm:.1f} nm"
    strain_text = f"microstrain: {microstrain:.1f}"

    return dict(
        x=x, y=y, sy=sy, y_fit=y_fit,
        grain_text=grain_text,
        strain_text=strain_text,
        label=label
    )


##############################################################################
# 1) Below we define the “instrument broadening” functions for each dataset.
##############################################################################

def y_calc_2021(df):
    """
    For 2021 data, you used:
        intermed = (Param_3_Value^2) - [(0.000811175888*(Center^2) +
                                         0.000022872696*Center +
                                         0.000308591696)^2]
        y = sqrt(intermed)*cos(Center)
    """
#NSLS-II for Ho2Ti2O7 0.000811175888x2 + 0.000022872696x + 0.000308591696
    # APS for HO2Ti2O7 -0.000291105820x2 + 0.000037751154x + 0.001301856675
    intermed = (
        (df['Param_3_Value']**2) -
        (
            (0.000811175888 * df['Center']**2) +
            ( 0.000022872696 * df['Center'])    +
            0.000308591696
        )**2
    )
    return np.sqrt(intermed) * np.cos(df['Center'])


def y_calc_2024(df):
#NSLS-II for 5 comp y = 0.00063190x2 + 0.00004343x + 0.00039145
#APS for 5 comp y = -0.00012982x2 - 0.00013165x + 0.00133723
    intermed = (
        df['Param_3_Value']**2 -
        ((0.00063190 * df['Center']**2)
         +(0.00004343 * df['Center'])
         + 0.00039145) ** 2
    )
    return np.sqrt(intermed) * np.cos(df['Center'])


##############################################################################
# 2) Provide file paths & process them
##############################################################################
#Pristine
file1 = "/Users/georgeadamson/Desktop/UTK python programs for reaserch/Williamson_Hall_peaks/Pristine GSI 2020/Ho_Ti/Ho2Ti2O7_pristine_20210617-121226_20e625_0001_mean_tth.peaks"
file2 = "/Users/georgeadamson/Desktop/UTK python programs for reaserch/Williamson_Hall_peaks/Pristine GSI 2020/HE_Ti/xrd_YbTbGdErDy2Ti2O7_HE_20231108-040640_7f99b9_primary-1_mean_tth.peaks"
#1E12
#file1 = '/Users/georgeadamson/Desktop/UTK python programs for reaserch/Williamson_Hall_peaks/Pristine GSI 2020/Ho_Ti/Ho2Ti2O7_1E12_001_mar2300_A0_bragg_peaks.peaks'
#file2 = '/Users/georgeadamson/Desktop/UTK python programs for reaserch/Williamson_Hall_peaks/Pristine GSI 2020/HE_Ti/High_Entropy_Titanate_1e12_001 Bragg Peaks.peaks'
d1 = process_peaks_file(file1, label="Ho2Ti2O7", y_calc_fn=y_calc_2021)
d2 = process_peaks_file(file2, label="(YbTbGdErDy)2Ti2O7", y_calc_fn=y_calc_2024)

##############################################################################
# 3) Plot all lines on the same graph
##############################################################################

plt.figure(figsize=(10, 6))
plt.tight_layout()

# Plot dataset #1
plt.errorbar(d1["x"], d1["y"]*1000, yerr=d1["sy"]*1000, fmt='o', label=f"{d1['label']} data")
plt.plot(d1["x"], d1["y_fit"]*1000, label=f"{d1['label']} fit")

# Plot dataset #2
plt.errorbar(d2["x"], d2["y"]*1000, yerr=d2["sy"]*1000, fmt='s', label=f"{d2['label']} data")
plt.plot(d2["x"], d2["y_fit"]*1000, label=f"{d2['label']} fit")

slope_e, intercept_e = np.polyfit(d2["x"], d2["y"], 1)
print("Excel‑style (unweighted) slope, intercept:", slope_e, intercept_e)
# Optionally annotate
plt.text(0.02, 0.77,
         f"{'Ho2Ti2O7'}\n{d1['grain_text']}\n{d1['strain_text']}",
         transform=plt.gca().transAxes, va='top',
         bbox=dict(boxstyle="round", facecolor="white"))
plt.text(0.02, 0.65,
         f"{'(YbTbGdErDy)2Ti2O7'}\n{d2['grain_text']}\n{d2['strain_text']}",
         transform=plt.gca().transAxes, va='top',
         bbox=dict(boxstyle="round", facecolor="white"))

plt.xlabel(r"$4\sin(\theta)$", fontsize=18)
plt.ylabel(r"$\beta\cos(\theta)$", fontsize=18)
plt.title('Pristine', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend()
plt.show()

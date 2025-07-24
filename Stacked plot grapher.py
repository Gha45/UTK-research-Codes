import os
import numpy as np
import matplotlib.pyplot as plt
import re
import pandas as pd
import sympy
# for exclusion range 1.3 to 4.7
SAVED_DIR_FILE = 'saved_directory.txt'
SAVED_EXCLUSION_RANGES_FILE = 'saved_exclusion_ranges.txt'


def get_plotting_range():
    """Prompt the user to input a desired x-range for plotting."""
    try:
        range_input = input(
            "Enter the desired x-range for plotting in the format 'start,end' (or press enter to use full range): ").strip()
        if not range_input:
            return None, None  # this means the full range will be used

        start, end = map(float, range_input.split(','))
        return start, end

    except ValueError:
        print("Invalid format. Using the full range.")
        return None, None


def normalize_data(data):
    min_val = np.min(data)
    max_val = np.max(data)

    return (data - min_val) / (max_val - min_val)


def get_directory_path():
    # Check if saved directory file exists
    if os.path.exists(SAVED_DIR_FILE):
        with open(SAVED_DIR_FILE, 'r') as file:
            saved_dir = file.readline().strip()

        # Ask the user if they want to use the saved directory
        choice = input(f"Do you want to use the saved directory '{saved_dir}'? (yes/no): ").strip().lower()
        if choice == 'yes':
            return saved_dir

    # If the user chooses not to use the saved directory or if no saved directory exists
    directory_path = input("Enter the directory path containing the data files: ")

    # Save the new directory for future reference
    with open(SAVED_DIR_FILE, 'w') as file:
        file.write(directory_path)

    return directory_path

# Use the function to get the directory path
directory_path = get_directory_path()


filenames = [f for f in os.listdir(directory_path) if f.endswith(('.xy', '.txt'))]
pristine_keywords = ['stock', 'unirr', 'original', 'clean','Pristine','pris']
def extract_numeric_value(filename):
    if "unirr" in filename or "stock" in filename:
        return (0, 0)

    match = re.search(r"(\d+(?:\.\d+)?)e(\d+)", filename)
    if match:
        base, exponent = match.groups()
        return (int(base), int(exponent))
    else:
        return (float('inf'), 0)


#sorted_filenames = sorted(filenames, key=extract_numeric_value)
sorted_filenames_and_values = sorted([(f, extract_numeric_value(f)) for f in filenames],
                                     key=lambda x: x[1][0] * 10**x[1][1])

# Ask the user for the exclusion ranges
def get_exclusion_ranges():
    # Check if saved exclusion ranges file exists
    if os.path.exists(SAVED_EXCLUSION_RANGES_FILE):
        with open(SAVED_EXCLUSION_RANGES_FILE, 'r') as file:
            saved_ranges = [tuple(map(float, line.strip().split(','))) for line in file.readlines()]
        choice = input(f"Do you want to use the saved exclusion ranges {saved_ranges}? (yes/no): ").strip().lower()
        if choice == 'yes':
            return saved_ranges

    # If the user chooses not to use the saved exclusion ranges or if no saved ranges exists
    exclude_ranges = []
    try:
        num_ranges = int(input("How many exclusion ranges do you want to input? (Enter 0 if none): "))
        if num_ranges > 0:
            for i in range(num_ranges):
                while True:
                    excl_range = input(f"Enter exclusion range {i + 1} in format 'start,end': ")
                    try:
                        start, end = map(float, excl_range.split(','))
                        exclude_ranges.append((start, end))
                        break
                    except ValueError:
                        print("Invalid format. Please enter the range as 'start,end'.")
    except ValueError:
        print("Invalid input. Assuming no exclusion ranges.")
        num_ranges = 0

    # Save the new exclusion ranges for future reference
    with open(SAVED_EXCLUSION_RANGES_FILE, 'w') as file:
        for r in exclude_ranges:
            file.write(f"{r[0]},{r[1]}\n")

    return exclude_ranges


cmap = plt.get_cmap("gist_heat")
exclude_ranges = get_exclusion_ranges()
plot_start, plot_end = get_plotting_range()
plot_title = input("Please enter the title for the plot: ")
offset = 1.0

# Constants for adjusting text label position
x_offset = 0.05 # Adjust this to set how far from the rightmost x-point the label should appear
#x_offset = 0.07
#y_position_factor = 1.1  # Adjust this to set how high the label should appear relative to the plot
# comment out if you want standard orientation
plt.figure(figsize=(8,10)) #8,14 8,10
for i, (f, (base, exponent)) in enumerate(sorted_filenames_and_values):
    full_path = os.path.join(directory_path, f)

    try:
        data = np.genfromtxt(full_path, encoding = 'latin 1')[:, :2]

        # First, filter the data based on the desired plotting range
        if plot_start is not None and plot_end is not None:
            data = data[(data[:, 0] >= plot_start) & (data[:, 0] <= plot_end)]

        # Then, create the mask for the exclusion ranges
        mask = np.ones(data[:, 0].shape, dtype=bool)
        for r in exclude_ranges:
            mask &= ~((data[:, 0] >= r[0]) & (data[:, 0] <= r[1]))

        # After filtering, apply the mask
        valid_x = data[:, 0][mask]
        valid_y = data[:, 1][mask]

        normalized_y = normalize_data(valid_y) + i*offset # Normalize data here
        average_y = np.mean(normalized_y)

        # Plot the data with adjusted line width and label font size
        color = cmap(i / (len(sorted_filenames_and_values) + 2))
        plt.plot(valid_x, normalized_y, label=f, linewidth=3.0,color=color)

        # Determine the rightmost x-value and the corresponding y-value
        max_x = valid_x[-1]
        corresponding_y = normalized_y[-1]
        if plot_start is None or plot_end is None:
            # Using the full range of data if plot_start or plot_end is None.
            plot_start = min(valid_x)
            plot_end = max(valid_x)

        if any(keyword in f for keyword in pristine_keywords):
            label_text = "Pristine"  # Set label_text to 'Pristine' if keyword is found
        else:
            # Otherwise, use the original method to set label_text
            label_text = f"{base}" + "E" + str(exponent)   # Original label_text format
        #plt.text(max_x + x_offset * (plot_end - plot_start), average_y, label_text,
          #       verticalalignment='center', fontsize= 20)
        plt.annotate(label_text, (max(valid_x)*0.8, i*offset+0.55), fontsize=30)
# 5.7 for Q data
    except IOError:
        print(f"Error reading file {f}. It might not exist or there's another issue.")
#plt.text(0.96, 0.98, '$(YbTbGdErDy)_2Ti_2O_7$', horizontalalignment='right',
#         verticalalignment='top', transform=plt.gca().transAxes, fontsize=20)
#plt.title(plot_title, fontsize = 24)
# $(YbTbGdErDy)_2Ti_2O_7$
#$(YbGd)_2Zr_2O_7$
#Ho$_2$Ti$_2$O$_7$
plt.title("SCP", fontsize = 30)

#plt.title("(YbTbGdErDy)$_2$Ti$_2$O$_7$", fontsize = 30)
#plt.xlabel("2"+r"${\Theta}$(Degrees)",fontsize= 24)
#plt.xlabel("Q ("+'$\AA\ ^{-1}$'+')',fontsize= 30)
plt.xlabel("Raman shift (cm$^-$$^1$)", fontsize = 30)
plt.xticks(fontsize = 24)
plt.yticks([])
plt.ylabel("Normalized Intensity (arb. units)",fontsize = 30)
plt.subplots_adjust( right=0.7,bottom = 0.15)
plt.tight_layout()
plt.savefig(plot_title+'.svg')
#plt.savefig("HE Zirc")
plt.show()

import numpy
import os
import math

# APS(in angstroms)
#  APS of 2025 :0.4239
#wavelength= 0.4133
wavelength = 0.4239
# NSLS-2
#wavelength = 0.1665

def modify_x_values(x):
    theta = x*(math.pi/180)
    return ((4* math.pi)/(wavelength))* math.sin(theta/2)

def process_files(source_folder, destination_folder):
    for filename in os.listdir(source_folder):
        if filename.endswith(".xy"):
            source_file = os.path.join(source_folder, filename)
            destination_file = os.path.join(destination_folder, filename)

            with open(source_file, 'r') as file:
                lines = file.readlines()

            modified_lines = []
            for line in lines:
                parts = line.split()
                if len(parts) == 2:
                    try:
                        x = float(parts[0])
                        y = parts[1]
                        x = modify_x_values(x)
                        modified_lines.append(f"{x} {y}\n")
                    except ValueError:
                        # Handle the case where conversion to float fails
                        modified_lines.append(line)

            with open(destination_file, 'w') as file:
                file.writelines(modified_lines)


source_folder = '/Users/georgeadamson/Desktop/UTK python programs for reaserch/Q data/Ho_APS2025'
destination_folder = '/Users/georgeadamson/Desktop/UTK python programs for reaserch/Q data/QHET'
process_files(source_folder, destination_folder)
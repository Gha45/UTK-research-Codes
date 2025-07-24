import os
import pandas as pd
import re
import math


def process_peak_files(input_folder, output_folder):
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Define a pattern to match "NUMBER +/- NUMBER"
    pattern = re.compile(r'([\d\.]+) \+/- ([\d\.e\+-]+)')

    # Define column names
    columns = ['PeakType', 'Center', 'Height', 'Area', 'FWHM'] + [f'Param_{i}{suffix}' for i in range(1, 5) for suffix
                                                                  in ('_Value', '_Error')]

    # Initialize an empty DataFrame to hold all results
    all_data = pd.DataFrame()

    # Traverse the directory
    for root, _, files in os.walk(input_folder):
        for file_name in files:
            if file_name.endswith('.peaks'):
                file_path = os.path.join(root, file_name)
                print("Processing:", file_name)
                # Read and process the file
                with open(file_path, 'r') as file:
                    lines = file.readlines()
                data_gaussian = []
                data_pseudovoigt = []
                for line in lines:
                    parts = line.strip().split()
                    if len(parts) > 5:
                        peak_type = parts[1]
                        center = parts[2]
                        height = parts[3]
                        area = parts[4]
                        fwhm = parts[5]
                        parameters = ' '.join(parts[6:])
                        matches = pattern.findall(parameters)
                        parsed_data = [peak_type, center, height, area, fwhm] + [m for match in matches for m in
                                                                                 match] + [None] * (8 - 2 * len(matches))

                        if "Gaussian" in line:
                            data_gaussian.append(parsed_data)
                        elif "PseudoVoigt" in line:
                            data_pseudovoigt.append(parsed_data)

                gaussian_df = pd.DataFrame(data_gaussian, columns=columns) if data_gaussian else pd.DataFrame(
                    {col: [0] * 1 for col in columns})
                pseudovoigt_df = pd.DataFrame(data_pseudovoigt, columns=columns) if data_pseudovoigt else pd.DataFrame(
                    {col: [0] * 1 for col in columns})
                # Remove Moly peaks if there
                pseudovoigt_df['Center'] = pd.to_numeric(pseudovoigt_df['Center'], errors='coerce')
                exclusion_ranges = [(10.65, 10.7), (15.0, 15.2)]
                for lower, upper in exclusion_ranges:
                    pseudovoigt_df = pseudovoigt_df[~pseudovoigt_df['Center'].between(lower, upper)]

                # Apply Mathematical equations to get amorph frac and error
                Amorph = gaussian_df['Area'].astype(float).sum()
                crystalline = pseudovoigt_df['Area'].astype(float).sum()
                Total = Amorph + crystalline
                amorph_frac = Amorph / Total
                amorph_error = math.sqrt(math.pow((gaussian_df['Param_1_Error'].astype(float).sum() / Total), 2) +
                                         math.pow(((Amorph * (math.sqrt(
                                             math.pow(gaussian_df['Param_1_Error'].astype(float).sum(), 2) + math.pow(
                                                 pseudovoigt_df['Param_1_Error'].astype(float).sum(), 2)))) / math.pow(
                                             Total, 2)), 2))

                final_df = pd.DataFrame({'File': file_name, 'Fraction': [amorph_frac], 'Error': [amorph_error]})

                # Append to the all_data DataFrame
                all_data = pd.concat([all_data, final_df], ignore_index=True)


    # Save the consolidated DataFrame to CSV
    output_csv_file = os.path.join(output_folder, "consolidated_peak_data.csv")
    all_data.to_csv(output_csv_file, index=False)
    print(f"All data saved to {output_csv_file}")

file_path = '/Users/georgeadamson/Desktop/Written_papers/Effects_of_SHI_in_CCO/GSAS-II Reitveld/Data/1-D_APS_2025/fityk_fit'
output_folder ='/Users/georgeadamson/Desktop/Written_papers/Effects_of_SHI_in_CCO/GSAS-II Reitveld/Data/APS_2025_fit'
process_peak_files(file_path, output_folder)

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
    columns = ['PeakType', 'Center', 'Height', 'Area', 'FWHM'] + [f'Param_{i}{suffix}' for i in range(1, 5) for suffix in ('_Value', '_Error')]

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
                        parsed_data = [peak_type, center, height, area, fwhm] + [m for match in matches for m in match] + [None] * (8 - 2 * len(matches))

                        if "Gaussian" in line:
                            data_gaussian.append(parsed_data)
                        elif "PseudoVoigt" in line:
                            data_pseudovoigt.append(parsed_data)

                # Create DataFrames
                gaussian_df = pd.DataFrame(data_gaussian, columns=columns) if data_gaussian else pd.DataFrame({col: [0] * 1 for col in columns})
                pseudovoigt_df = pd.DataFrame(data_pseudovoigt, columns=columns) if data_pseudovoigt else pd.DataFrame({col: [0] * 1 for col in columns})

                # Filter PseudoVoigt DataFrame
                pseudovoigt_df['Center'] = pd.to_numeric(pseudovoigt_df['Center'], errors='coerce')
                exclusion_ranges = [(10.64, 10.7), (15.0, 15.2)]
                for lower, upper in exclusion_ranges:
                    pseudovoigt_df = pseudovoigt_df[~pseudovoigt_df['Center'].between(lower, upper)]

                # Calculate amorphous fraction and error
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

    # Calculating combined error
    def extract_identifier(file_name):
        match = re.search(r'(\d+E\d+|pris)', file_name)
        return match.group(0) if match else None

    all_data['Identifier'] = all_data['File'].apply(extract_identifier)
    grouped_data = all_data.groupby('Identifier').agg({
        'Fraction': 'mean',
        'Error': lambda x: math.sqrt(sum(x**2))
    }).reset_index()
    grouped_data.columns = ['Identifier', 'Average Fraction', 'Combined Error']

    # Save the consolidated and grouped data to CSV
    output_csv_file = os.path.join(output_folder, "peak_data.csv")
    all_data.to_csv(output_csv_file, index=False)
    print(f"All data saved to {output_csv_file}")

    grouped_csv_file = os.path.join(output_folder, "grouped_peak_data.csv")
    grouped_data.to_csv(grouped_csv_file, index=False)
    print(f"Grouped data saved to {grouped_csv_file}")


#input file path for the folder full of the peak files
#file_path = '/Users/georgeadamson/Desktop/Data/Cades data analysis/Archive/HZTO_Peaks'
# output file  for the csv files
#output_folder = '/Users/georgeadamson/Desktop/Data/Cades data analysis/Archive/Testing_HZTO'

file_path = '/Users/georgeadamson/Desktop/Data/1-D_APS_2025/HTO-test'
output_folder ='/Users/georgeadamson/Desktop/Data/1-D_APS_2025/HTO-testCase'
process_peak_files(file_path, output_folder)

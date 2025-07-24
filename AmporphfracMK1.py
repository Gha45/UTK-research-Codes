import os
import pandas as pd
import re
import math
def process_peak_files(input_folder, output_folder):
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Define a pattern to match "NUMBER +/- NUMBER"
    pattern = re.compile(r'([\d\.]+) \+/- ([\d\.e\+-]+)')

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
                        peak_type = parts[1]  # Assumes the first part is a marker like '%_1'
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


                # Define column names
                columns = ['PeakType', 'Center', 'Height', 'Area', 'FWHM'] + [f'Param_{i}{suffix}' for i in range(1, 5) for suffix in ('_Value', '_Error')]

                # Create DataFrames
                gaussian_df = pd.DataFrame(data_gaussian, columns=columns) if data_gaussian else pd.DataFrame(
                    {col: [0] * 1 for col in columns})
                pseudovoigt_df = pd.DataFrame(data_pseudovoigt, columns=columns) if data_pseudovoigt else pd.DataFrame(
                    {col: [0] * 1 for col in columns})

                # Apply mathematical transformations (example: adding a new column for normalized height
                Amorph = (gaussian_df['Area'].astype(float).sum())
                crystalie = (pseudovoigt_df['Area'].astype(float).sum())

                Total = Amorph + crystalie
                amorph_frac = (Amorph)/(Total)

                amorph_error = math.sqrt(math.pow((gaussian_df['Param_1_Error'].astype(float).sum() / Total),2) +
                                          math.pow(((Amorph *(math.sqrt(math.pow(gaussian_df['Param_1_Error'].astype(float).sum(),2)+ math.pow(pseudovoigt_df['Param_1_Error'].astype(float).sum(),2))))/ math.pow(Total,2)),2))

                final_df = pd.DataFrame({'Fraction': [amorph_frac], 'Error': [amorph_error]})
                print(final_df)
                # Save each processed file to a separate output file
                output_file = os.path.join(output_folder, f"{os.path.splitext(file_name)[0]}_processed.xlsx")
                final_df.to_excel(output_file, index=False)

                print(f"Data saved to {output_file}")

file_path = '/Users/georgeadamson/Desktop/Data/Cades data analysis/Archive/HZTO_Peaks'
output_folder ='/Users/georgeadamson/Desktop/Data/Cades data analysis/Archive/Testing_HZTO'
process_peak_files(file_path, output_folder)

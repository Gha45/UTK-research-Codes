import os
import pandas as pd
import re
import math

def process_peak_files(input_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)

    pattern = re.compile(r'([\d\.]+) \+/- ([\d\.e\+-]+)')
    columns = ['PeakType', 'Center', 'Height', 'Area', 'FWHM'] + [f'Param_{i}{suffix}' for i in range(1, 5) for suffix in ('_Value', '_Error')]

    all_data = pd.DataFrame()

    for root, _, files in os.walk(input_folder):
        for file_name in files:
            if file_name.endswith('.peaks'):
                file_path = os.path.join(root, file_name)
                print("Processing:", file_name)

                with open(file_path, 'r') as file:
                    lines = file.readlines()

                data_gaussian, data_pseudovoigt = [], []

                for line in lines:
                    parts = line.strip().split()
                    if len(parts) > 5:
                        peak_type, center, height, area, fwhm = parts[1:6]
                        parameters = ' '.join(parts[6:])
                        matches = pattern.findall(parameters)
                        parsed_data = [peak_type, center, height, area, fwhm] + [m for match in matches for m in match] + [None] * (8 - 2 * len(matches))

                        if "Gaussian" in line:
                            data_gaussian.append(parsed_data)
                        elif "PseudoVoigt" in line:
                            data_pseudovoigt.append(parsed_data)

                gaussian_df = pd.DataFrame(data_gaussian, columns=columns) if data_gaussian else pd.DataFrame({col: [0] * 1 for col in columns})
                pseudovoigt_df = pd.DataFrame(data_pseudovoigt, columns=columns) if data_pseudovoigt else pd.DataFrame({col: [0] * 1 for col in columns})

                pseudovoigt_df['Center'] = pd.to_numeric(pseudovoigt_df['Center'], errors='coerce')
                exclusion_ranges = []  # [(10.64, 10.7), (15.0, 15.2)]
                for lower, upper in exclusion_ranges:
                    pseudovoigt_df = pseudovoigt_df[~pseudovoigt_df['Center'].between(lower, upper)]

                Amorph = gaussian_df['Area'].astype(float).sum()
                crystalline = pseudovoigt_df['Area'].astype(float).sum()
                Total = Amorph + crystalline

                # Calculate amorphous fraction
                amorph_frac = Amorph / Total if Total > 0 else 0

                # Handle missing errors
                gaussian_df['Param_1_Error'] = gaussian_df['Param_1_Error'].astype(float).fillna(0)
                pseudovoigt_df['Param_1_Error'] = pseudovoigt_df['Param_1_Error'].astype(float).fillna(0)

                # Error propagation
                if Total > 0:
                    if Amorph == 0:
                        # No amorphous peaks, fraction = 0, error depends on crystalline errors
                        amorph_error = math.sqrt(pseudovoigt_df['Param_1_Error'].pow(2).sum()) / Total
                    elif crystalline == 0:
                        # No crystalline peaks, fraction = 1, error depends on amorphous errors
                        amorph_error = math.sqrt(gaussian_df['Param_1_Error'].pow(2).sum()) / Total
                    else:
                        # Both amorphous and crystalline peaks exist
                        amorph_error = math.sqrt(
                            (crystalline ** 2 * gaussian_df['Param_1_Error'].pow(2).sum() +
                             Amorph ** 2 * pseudovoigt_df['Param_1_Error'].pow(2).sum()) / Total ** 4
                        )
                else:
                    # No peaks at all
                    amorph_error = 0

                final_df = pd.DataFrame({'File': [file_name], 'Fraction': [amorph_frac], 'Error': [amorph_error]})
                all_data = pd.concat([all_data, final_df], ignore_index=True)

    def extract_identifier(file_name):
        match = re.search(r'(\d+[eE]\d+)|pris', file_name)
        if match:
            identifier = match.group(0)
            return identifier.upper() if identifier != 'pris' else identifier
        return None

    def identifier_to_sort_key(identifier):
        if identifier == 'pris':
            return float('inf')  # Place 'pris' at the end
        try:
            return float(identifier.replace('E', 'e'))  # Convert to float (e.g., "1E23" -> 1e23)
        except ValueError:
            return float('inf')  # Fallback for invalid identifiers

    all_data['Identifier'] = all_data['File'].apply(extract_identifier)
    grouped_data = all_data.groupby('Identifier').agg({
        'Fraction': 'mean',
        'Error': lambda x: math.sqrt(sum(x ** 2))
    }).reset_index()
    grouped_data.columns = ['Identifier', 'Average Fraction', 'Combined Error']

    # Sort grouped_data by numerical value of Identifier
    grouped_data['SortKey'] = grouped_data['Identifier'].apply(identifier_to_sort_key)
    grouped_data = grouped_data.sort_values('SortKey').drop('SortKey', axis=1)

    output_csv_file = os.path.join(output_folder, "peak_data.csv")
    all_data.to_csv(output_csv_file, index=False)
    print(f"All data saved to {output_csv_file}")

    grouped_csv_file = os.path.join(output_folder, "grouped_peak_data.csv")
    grouped_data.to_csv(grouped_csv_file, index=False)
    print(f"Grouped data saved to {grouped_csv_file}")

# Example usage
file_path = '/Users/georgeadamson/Desktop/Amorphous Fraction anaylsis/Fityk fitting/Complete Ho2Ti2O7'
output_folder = '/Users/georgeadamson/Desktop/Amorphous Fraction anaylsis/Amorphrac3 Python Script output/Ho2Ti2O7_Output3'
process_peak_files(file_path, output_folder)
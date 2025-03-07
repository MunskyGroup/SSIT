import pandas as pd
import numpy as np

# Load the CSV file
file_path = 'pdoCalibrationData_NeuertLab_240925.csv'
data = pd.read_csv(file_path)

# Filter the data based on Condition = '2M_NaCl_Step'
filtered_data = data.loc[data['Condition'].str.contains('2M_NaCl_Step', case=False, na=False), 
                         ['Cell_id', 'Condition', 'Replica', 'Conc_mM', 'Time_index_min', 
                          'RNA_STL1_cyto_TS3Full', 'RNA_STL1_nuc_TS3Full', 'Cyto_STL1_avg_int_TS3Full', 'Nuc_STL1_avg_int_TS3Full']]

# Sum 'RNA_STL1_cyto_TS3Full' and 'RNA_STL1_nuc_TS3Full' into new column 'RNA_STL1_total_TS3Full'
filtered_data['RNA_STL1_total_TS3Full'] = filtered_data['RNA_STL1_cyto_TS3Full'] + filtered_data['RNA_STL1_nuc_TS3Full']

# Sum 'Cyto_STL1_avg_int_TS3Full' and 'Nuc_STL1_avg_int_TS3Full' into new column 'STL1_avg_int_TS3Full' and transform it into an integer for SSIT
#filtered_data['STL1_avg_int_TS3Full'] = (filtered_data['Cyto_STL1_avg_int_TS3Full'] + filtered_data['Nuc_STL1_avg_int_TS3Full']).astype(int)
# Fill NaN or inf values with 0, then sum and convert to integer
filtered_data['STL1_avg_int_TS3Full'] = (
    filtered_data['Cyto_STL1_avg_int_TS3Full'].replace([np.inf, -np.inf], np.nan).fillna(0) +
    filtered_data['Nuc_STL1_avg_int_TS3Full'].replace([np.inf, -np.inf], np.nan).fillna(0)
).astype(int)


# Replace 'Time_index_min' with 'time' to be read by SSIT.m
filtered_data.rename(columns={'Time_index_min': 'time'}, inplace=True)

# Define the output file path
output_file = 'filtered_data_2M_NaCl_Step.csv'

# Save the filtered data to a new CSV file
filtered_data.to_csv(output_file, index=False)

output_file

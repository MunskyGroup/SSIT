import pandas as pd

# Load the CSV file
file_path = 'pdoCalibrationData_NeuertLab_240925.csv'
data = pd.read_csv(file_path)

# Filter the data based on Condition = '2M_NaCl_Step'
# Adjusting the filter to be case-insensitive and checking for partial matches
filtered_data = data.loc[data['Condition'].str.contains('2M_NaCl_Step', case=False, na=False), 
                         ['Cell_id', 'Condition', 'Replica', 'Conc_mM', 'time', 'RNA_STL1_cyto_TS3Full']]
                         
# Define the output file path
output_file_path = 'filtered_data_2M_NaCl_Step.csv'

# Save the filtered data to a new CSV file
filtered_data.to_csv(output_file_path, index=False)

output_file_path

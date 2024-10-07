import pandas as pd

# Load the CSV file
file_path = 'pdoCalibrationData_NeuertLab_240925.csv'
data = pd.read_csv(file_path)

# Filter the data based on Condition = '2M_NaCl_Step'
filtered_data = data.loc[data['Condition'].str.contains('2M_NaCl_Step', case=False, na=False), 
                         ['Cell_id', 'Condition', 'Replica', 'Conc_mM', 'Time_index_min', 'RNA_STL1_cyto_TS3Full']]

filtered_data.rename(columns={'Time_index_min': 'time'}, inplace=True)


# Define the output file path
output_file = 'filtered_data_2M_NaCl_Step.csv'


# Save the filtered data to a new CSV file
filtered_data.to_csv(output_file, index=False)

output_file

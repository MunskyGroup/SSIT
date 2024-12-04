import pandas as pd

# Load the CSV files
file1_path = 'dataframe_A549_DUSP1_100nM_10min_062723_spot_counts.csv'
file2_path = 'dataframe_A549_DUSP1_100nM_10min_062723_spot_counts_z_5_9.csv'
file3_path = 'dataframe_A549_DUSP1_100nM_10min_062723_spot_counts_z_7.csv'

df1 = pd.read_csv(file1_path)
df2 = pd.read_csv(file2_path)
df3 = pd.read_csv(file3_path)

# Rename the "Spot_Count" column in the second dataframe to "Z_5_9"
df2.rename(columns={'Spot_Count': 'Z_5_9'}, inplace=True)

# Rename the "Spot_Count" column in the third dataframe to "Z_7"
df3.rename(columns={'Spot_Count': 'Z_7'}, inplace=True)

# Merge the three dataframes on the "Image_Cell_ID" column
df4 = pd.merge(df1, df2, on='Image_Cell_ID', how='inner')
merged_df = pd.merge(df4, df3, on='Image_Cell_ID', how='inner')

# Save the merged dataframe to a new CSV file
merged_file_path = 'merged_dataframe_A549_DUSP1_100nM_10min_062723.csv'
merged_df.to_csv(merged_file_path, index=False)

merged_file_path

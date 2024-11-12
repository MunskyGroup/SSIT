# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt

# File paths
file_paths = [
    'dataframe_A549_DUSP1_100nM_10min_062723.csv',
    'dataframe_A549_DUSP1_100nM_20min_062723.csv',
    'dataframe_A549_DUSP1_100nM_30min_062723.csv',
    'dataframe_A549_DUSP1_100nM_40min_062723.csv',
]

# Corresponding time points
time_points = ['10min', '20min', '30min', '40min']

# Initialize dictionaries to store spot counts
total_spot_counts_by_time = {}
spot_counts_z_0_3_by_time = {}
spot_counts_z_4_7_by_time = {}
spot_counts_z_8_11_by_time = {}
spot_counts_z_12_15_by_time = {}
spot_counts_z_16_19_by_time = {}
spot_counts_z_20_23_by_time = {}
spot_counts_z_24_27_by_time = {}

# Loop through each file and calculate the spot counts
for file, time_point in zip(file_paths, time_points):
    df = pd.read_csv(file)
    
    # Count unique spot_id values per image_id and cell_id for all z-values
    total_spot_count = df.groupby(['image_id', 'cell_id'])['spot_id'].nunique().sum()
    
    # Filter for z-slices
    filtered_df_z_0_3 = df[(df['z'] >= 0) & (df['z'] <= 3)]
    filtered_df_z_4_7 = df[(df['z'] >= 4) & (df['z'] <= 7)]
    filtered_df_z_8_11 = df[(df['z'] >= 8) & (df['z'] <= 11)]
    filtered_df_z_12_15 = df[(df['z'] >= 12) & (df['z'] <= 15)]
    filtered_df_z_16_19 = df[(df['z'] >= 16) & (df['z'] <= 19)]
    filtered_df_z_20_23 = df[(df['z'] >= 20) & (df['z'] <= 23)]
    filtered_df_z_24_27 = df[(df['z'] >= 24) & (df['z'] <= 27)]
    
    # Count unique spot_id values per image_id and cell_id for z=0-3
    spot_count_z_0_3 = filtered_df_z_0_3.groupby(['image_id', 'cell_id'])['spot_id'].nunique().sum()
    
    # Count unique spot_id values per image_id and cell_id for z=4-7
    spot_count_z_4_7 = filtered_df_z_4_7.groupby(['image_id', 'cell_id'])['spot_id'].nunique().sum()
    
    # Count unique spot_id values per image_id and cell_id for z=8-11
    spot_count_z_8_11 = filtered_df_z_8_11.groupby(['image_id', 'cell_id'])['spot_id'].nunique().sum()
    
    # Count unique spot_id values per image_id and cell_id for z=12-15
    spot_count_z_12_15 = filtered_df_z_12_15.groupby(['image_id', 'cell_id'])['spot_id'].nunique().sum()
    
    # Count unique spot_id values per image_id and cell_id for z=16-19
    spot_count_z_16_19 = filtered_df_z_16_19.groupby(['image_id', 'cell_id'])['spot_id'].nunique().sum()
    
    # Count unique spot_id values per image_id and cell_id for z=20-23
    spot_count_z_20_23 = filtered_df_z_20_23.groupby(['image_id', 'cell_id'])['spot_id'].nunique().sum()
    
     # Count unique spot_id values per image_id and cell_id for z=24-27
    spot_count_z_24_27 = filtered_df_z_24_27.groupby(['image_id', 'cell_id'])['spot_id'].nunique().sum()
    
    # Save the total spot counts for this time point
    total_spot_counts_by_time[time_point] = total_spot_count
    spot_counts_z_0_3_by_time[time_point] = spot_count_z_0_3
    spot_counts_z_4_7_by_time[time_point] = spot_count_z_4_7
    spot_counts_z_8_11_by_time[time_point] = spot_count_z_8_11
    spot_counts_z_12_15_by_time[time_point] = spot_count_z_12_15
    spot_counts_z_16_19_by_time[time_point] = spot_count_z_16_19
    spot_counts_z_20_23_by_time[time_point] = spot_count_z_20_23
    spot_counts_z_24_27_by_time[time_point] = spot_count_z_24_27

# Prepare data for plotting
times = list(total_spot_counts_by_time.keys())
total_spot_counts = list(total_spot_counts_by_time.values())
spot_counts_z_0_3 = list(spot_counts_z_0_3_by_time.values())
spot_counts_z_4_7 = list(spot_counts_z_4_7_by_time.values())
spot_counts_z_8_11 = list(spot_counts_z_8_11_by_time.values())
spot_counts_z_12_15 = list(spot_counts_z_12_15_by_time.values())
spot_counts_z_16_19 = list(spot_counts_z_16_19_by_time.values())
spot_counts_z_20_23 = list(spot_counts_z_20_23_by_time.values())
spot_counts_z_24_27 = list(spot_counts_z_24_27_by_time.values())

# Create the plot
plt.figure(figsize=(10,6))

# Plot total spot counts across all z-values
plt.plot(times, total_spot_counts, marker='o', linestyle='-', color='b', markersize=8, label='Total Spot Counts')

# Plot z=24-27 spot counts
plt.plot(times, spot_counts_z_24_27, marker='o', linestyle='-', color='c', markersize=8, label='Spot Counts (z=24-27)')

# Plot z=20-23 spot counts
plt.plot(times, spot_counts_z_20_23, marker='o', linestyle='-', color='c', markersize=8, label='Spot Counts (z=20-23)')

# Plot z=16-19 spot counts
plt.plot(times, spot_counts_z_16_19, marker='o', linestyle='-', color='c', markersize=8, label='Spot Counts (z=16-19)')

# Plot z=12-15 spot counts
plt.plot(times, spot_counts_z_12_15, marker='o', linestyle='-', color='c', markersize=8, label='Spot Counts (z=12-15)')

# Plot z=8-11 spot counts
plt.plot(times, spot_counts_z_8_11, marker='o', linestyle='-', color='c', markersize=8, label='Spot Counts (z=8-11)')

# Plot z=4-7 spot counts
plt.plot(times, spot_counts_z_4_7, marker='o', linestyle='-', color='r', markersize=8, label='Spot Counts (z=4-7)')

# Plot z=0-3 spot counts
plt.plot(times, spot_counts_z_0_3, marker='o', linestyle='-', color='r', markersize=8, label='Spot Counts (z=0-3)')


# Add labels and title
plt.title('Spot Counts Through Time (All z-stacks and z-slice subsets) 100nM Dex')
plt.xlabel('Time (min)')
plt.ylabel('Spot Count')
plt.grid(True)

# Adjust the legend placement to move it slightly up from the lower right corner
plt.legend(loc='upper left', bbox_to_anchor=(0, 1))

# Display the plot
plt.show()
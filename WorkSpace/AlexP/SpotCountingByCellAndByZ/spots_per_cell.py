import pandas as pd

# List of file names
files = [
    'dataframe_A549_DUSP1_100nM_10min_062723.csv',
    'dataframe_A549_DUSP1_100nM_20min_062723.csv',
    'dataframe_A549_DUSP1_100nM_30min_062723.csv',
    'dataframe_A549_DUSP1_100nM_40min_062723.csv'
]

# Initialize a dictionary to store the total "spot_count" for each file
spot_count_per_file = {}

# Loop through each file
for file in files:
    # Load the CSV file
    df = pd.read_csv(file)

    # Initialize a dictionary to store the spot count for each "cell_id" in the file
    file_spot_counts = {}

    # Loop through each "image_id"
    for image_id in df['image_id'].unique():
        image_df = df[df['image_id'] == image_id]

        # Loop through each "cell_id" in the current "image_id"
        for cell_id in image_df['cell_id'].unique():
            cell_df = image_df[image_df['cell_id'] == cell_id]

            # Count the number of unique "spot_id" values per "z" value for this cell_id
            spot_count = cell_df.groupby('z')['spot_id'].nunique().sum()

            # Save the total "spot_count" for this "cell_id"
            file_spot_counts[(image_id, cell_id)] = spot_count

    # Store the result for the current file
    spot_count_per_file[file] = file_spot_counts

# Save the spot counts for each file into separate CSVs
for file, spot_counts in spot_count_per_file.items():
    # Convert the dictionary to a DataFrame
    spot_counts_df = pd.DataFrame(list(spot_counts.items()), columns=['Image_Cell_ID', 'Spot_Count'])
    
    # Save the DataFrame to a new CSV file
    output_filename = file.replace('.csv', '_spot_counts.csv')
    spot_counts_df.to_csv(output_filename, index=False)

print("Spot count files saved successfully.")
import pandas as pd

# Load the CSV file
file_path = "GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03.csv"
df = pd.read_csv(file_path)

# Check the first few rows to understand the structure
df.head()

# Define the required columns
required_columns = ["nucGRint", "cytoGRint", "time", "dex_conc", "unique_cell_id", "normGRcyt", "normGRnuc"]

# Filter the DataFrame to only include the required columns
df_filtered = df[required_columns]

# Define the output file paths
output_files = {
    1: "GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03_dex_conc_1.csv",
    10: "GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03_dex_conc_10.csv",
    100: "GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03_dex_conc_100.csv"
}

# Filter and save the data for each Dex_Conc value
for conc, path in output_files.items():
    df_filtered_conc = df_filtered[df_filtered["dex_conc"] == conc]
    df_filtered_conc.to_csv(path, index=False)

# Provide confirmation
output_files
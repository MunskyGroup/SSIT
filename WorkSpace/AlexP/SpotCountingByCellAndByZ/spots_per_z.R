# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load the dataset (adjust the path as needed)
df <- read.csv("dataframe_A549_DUSP1_100nM_10min_062723.csv")

# Create a line plot for each image_id and cell_id
unique_image_ids <- unique(df$image_id)

for (image_id in unique_image_ids) {
  # Filter the data for the current image_id
  filtered_data_image <- df %>%
    filter(image_id == !!image_id)
  
  # Group by 'cell_id' and 'z', and count the unique 'spot_id' values for each combination
  unique_spots_per_z_cell <- filtered_data_image %>%
    group_by(cell_id, z) %>%
    summarize(unique_spot_count = n_distinct(spot_id))
  
  # Create the plot
  p <- ggplot(unique_spots_per_z_cell, aes(x = z, y = unique_spot_count, group = cell_id)) +
    geom_line(alpha = 0.5) +
    geom_point(alpha = 0.5) +
    labs(title = paste("Number of Unique spot_id Values for Each z by cell_id in image_id =", image_id),
         x = "z Values", 
         y = "Unique Spot Count") +
    theme_minimal()
  
  # Display the plot
  print(p)
}

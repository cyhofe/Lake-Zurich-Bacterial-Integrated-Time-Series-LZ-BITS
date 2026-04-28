# Clear R's environment
rm(list = ls())

# Load necessary libraries
library(utils)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)

# Set working directory
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Gene_Data/Raw_Data/KEGG_Categories/")

# Define the base folder containing KEGG category files
base_folder <- "KEGG_Categories/"

# List all major categories (top-level folders)
categories <- list.dirs(base_folder, recursive = FALSE)

# Initialize an empty dataframe with correct column names
kegg_df <- data.frame(KEGG_ID = character(), FunctionalCategory = character(), FunctionalSubcategory = character(), stringsAsFactors = FALSE)

# Loop through each category folder
for (category_path in categories) {
  
  # Extract category name from folder path
  category_name <- basename(category_path)
  
  # Handle Transporters category separately
  if (category_name == "Transporters") {
    
    # Get all second-level subdirectories (e.g., ABC, MFS, PTS)
    transporter_subcategories <- list.dirs(category_path, recursive = FALSE)
    
    for (subfolder in transporter_subcategories) {
      subcategory_name <- basename(subfolder)
      
      # List all .txt files inside this transport subcategory
      subcategory_files <- list.files(subfolder, full.names = TRUE, pattern = "\\.txt$")
      
      for (file_path in subcategory_files) {
        
        # Extract sub-subcategory name (file name without .txt extension)
        sub_subcategory_name <- basename(file_path) %>% str_remove("\\.txt$")
        
        # Read the KO entries (Ensures only KEGG_ID is extracted)
        ko_data <- read_delim(file_path, delim = "\t", col_names = c("KEGG_ID", "Ignore"), col_types = "cc") %>%
          select(KEGG_ID)  # Keep only KEGG_ID
        
        # Assign both category and subcategory properly
        ko_data <- ko_data %>%
          mutate(
            FunctionalCategory = category_name,
            FunctionalSubcategory = paste(subcategory_name, sub_subcategory_name, sep = " - ")  # Combine both levels
          )
        
        # Append to main dataframe
        kegg_df <- bind_rows(kegg_df, ko_data)
      }
    }
    
  } else {
    
    # Process other categories normally
    subcategory_files <- list.files(category_path, full.names = TRUE, pattern = "\\.txt$")
    
    for (file_path in subcategory_files) {
      
      # Extract subcategory name from filename
      subcategory_name <- basename(file_path) %>% str_remove("\\.txt$")
      
      # Read the KO entries (Ensures only KEGG_ID is extracted)
      ko_data <- read_delim(file_path, delim = "\t", col_names = c("KEGG_ID", "Ignore"), col_types = "cc") %>%
        select(KEGG_ID)  # Keep only KEGG_ID
      
      # Assign category and subcategory
      ko_data <- ko_data %>%
        mutate(
          FunctionalCategory = category_name,
          FunctionalSubcategory = subcategory_name
        )
      
      # Append to main dataframe
      kegg_df <- bind_rows(kegg_df, ko_data)
    }
  }
}

# Ensure final output contains only the required columns
kegg_df <- kegg_df %>% select(KEGG_ID, FunctionalCategory, FunctionalSubcategory)

# Save the final compiled dataset
output_file <- file.path("KEGG_Categories_Combined_LongFormat.csv")
write_csv(kegg_df, output_file)

# Print confirmation
print("✅ KEGG functional categories successfully compiled!")
print(head(kegg_df))  # Display first few rows

# unique ids and ids that are in multiple cateogries
kegg_df %>%
  summarise(
    Total_KEGG_IDs = n_distinct(KEGG_ID),  # Count unique KEGG IDs
    Duplicated_KEGG_IDs = sum(duplicated(KEGG_ID))  # Count duplicated KEGG IDs
  )


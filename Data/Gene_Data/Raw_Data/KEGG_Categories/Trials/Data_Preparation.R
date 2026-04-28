#clear R's brain
rm(list = ls())

#libraries
library(utils)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)

#set working directory
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Gene_Data/Raw_Data/KEGG_Categories/")

# brite categories #####################################################################################################################################

# Directory containing the files
brite_dir <- "kegg_brite_categories/"  # Replace with your actual directory path

# Get all file paths in the directory
brite_paths <- list.files(brite_dir, pattern = "\\.txt$", full.names = TRUE)

# Read files into a named list of dataframes
brite_list <- lapply(brite_paths, function(file) {
  read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)  # Adjust as needed
})

# Name each dataframe in the list by its filename (without extension)
names(brite_list) <- basename(brite_paths) %>% 
  tools::file_path_sans_ext()

# Check the resulting list
print(names(brite_list))  # Prints the names of the dataframes

# Combine all KEGG IDs into one unique list
all_kegg_ids <- unique(unlist(lapply(brite_list, function(df) df$V1)))

# Initialize a dataframe with all unique KEGG IDs
brite_categories_dat <- data.frame(KEGG_ID = all_kegg_ids, stringsAsFactors = FALSE)

# Loop through each dataframe in brite_list and add a column to brite_categories_dat
for (name in names(brite_list)) {
  # Check if each KEGG_ID is in the current dataframe
  brite_categories_dat[[name]] <- brite_categories_dat$KEGG_ID %in% brite_list[[name]]$V1
}

# Check the resulting dataframe
head(brite_categories_dat)

#export data
write.csv(brite_categories_dat, file = "Brite_KEGG_Categories.csv", row.names = FALSE)

# transporter categories ##############################################################################################################################

# Directory containing the files
transporter_dir <- "kegg_transporter_categories/"  # Replace with your actual directory path

# Get all subdirectories (folders)
folders <- list.dirs(transporter_dir, recursive = FALSE)

# Initialize an empty list to store all folders and their files
transporter_category_list <- list()

# Loop through each folder
for (folder in folders) {
  # Get the folder name
  folder_name <- basename(folder)
  
  # Get all .txt files in the folder
  files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
  
  # Read files into a named list of dataframes
  transporter_list <- lapply(files, function(file) {
    read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)  # Adjust as needed
  })
  
  # Name the dataframes in the list after the file names (without extension)
  names(transporter_list) <- tools::file_path_sans_ext(basename(files))
  
  # Add the list of dataframes to the transporter_category_list, named after the folder
  transporter_category_list[[folder_name]] <- transporter_list
}

# Check the resulting structure
str(transporter_category_list)

# Initialize an empty list to collect combined data
combined_data <- list()

# Loop through transporter_category_list
for (folder_name in names(transporter_category_list)) {
  folder_data <- transporter_category_list[[folder_name]]
  for (file_name in names(folder_data)) {
    df <- folder_data[[file_name]]
    df <- df %>%
      mutate(Folder = folder_name, File = file_name)  # Add folder and file info
    combined_data[[paste(folder_name, file_name, sep = "_")]] <- df
  }
}

# Combine all dataframes into one
final_combined_df <- bind_rows(combined_data, .id = "Category")

# Validate the uniqueness of KEGG_ID within each category
validation <- final_combined_df %>%
  group_by(V1, Folder, File) %>%
  summarize(count = n(), .groups = "drop") %>%
  filter(count > 1)

if (nrow(validation) > 0) {
  warning("Duplicate KEGG_IDs detected within categories. Investigate the following:")
  print(validation)
} else {
  message("No duplicate KEGG_IDs within categories.")
}

# Clean the combined data
transporter_df <- final_combined_df %>%
  rename("KEGG_ID" = "V1",
         "Transporter_Category" = "Folder",
         "Transporter_Type" = "File") %>%
  mutate(Transporter = paste(Transporter_Category, Transporter_Type, sep = "_")) %>% 
  select(KEGG_ID, Transporter) %>% distinct(KEGG_ID, Transporter, .keep_all = TRUE)

# Create wide format with TRUE/FALSE
wide_transporter_df <- transporter_df %>%
  mutate(Value = TRUE) %>%  # Add a column with TRUE
  pivot_wider(
    names_from = Transporter,  # Use column name directly
    values_from = Value,
    values_fill = FALSE  # Fill missing values with FALSE
  )

# Check the resulting dataframe
head(wide_transporter_df)

# Identify rows with multiple TRUEs
rows_with_multiple_trues <- wide_transporter_df %>%
  mutate(Multiple_TRUEs = rowSums(select(., -KEGG_ID)) > 1) %>%  # Calculate number of TRUEs per row
  filter(Multiple_TRUEs)  # Keep only rows with multiple TRUEs

# View the resulting rows
print(rows_with_multiple_trues)

#export data
write.csv(wide_transporter_df, file = "Transporter_KEGG_Categories.csv", row.names = FALSE)



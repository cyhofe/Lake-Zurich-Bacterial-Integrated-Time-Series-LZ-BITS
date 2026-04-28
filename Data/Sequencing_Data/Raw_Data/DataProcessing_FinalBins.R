# Clear R's brain
rm(list = ls())

# Libraries
library(utils)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)

# Set working directory
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Sequencing_Data/")

# Paths to directories
betula_path <- "Raw_Data/Betula_Workflow_Statistics/4_Binning/Final_Bins/"
ginkgo_path <- "Raw_Data/Ginkgo_Workflow_Statistics/4_Binning/Final_Bins/"
quercus_path <- "Raw_Data/Quercus_Workflow_Statistics/4_Binning/Final_Bins/"
salix_path <- "Raw_Data/Salix_Workflow_Statistics/4_Binning/Final_Bins/"

# read seqstat data ########################################################################################

# Function to fetch and read .txt files into a list of data frames
read_txt_files <- function(directory) {
  # Get all .txt file paths in the directory
  file_paths <- list.files(directory, pattern = "\\.txt$", full.names = TRUE)
  
  # Read each .txt file into a data frame
  data_list <- lapply(file_paths, function(file) {
    read_delim(file, delim = " ") # Force all columns to character
  })
  
  names(data_list) <- basename(file_paths) # Use file names as list names
  
  return(data_list)
}

# Load data from each directory
betula_seqstat <- read_txt_files(betula_path)
ginkgo_seqstat <- read_txt_files(ginkgo_path)
quercus_seqstat <- read_txt_files(quercus_path)
salix_seqstat <- read_txt_files(salix_path)

# Combine all lists into a single list of data frames
seqstat_data <- list(
  Betula = betula_seqstat,
  Ginkgo = ginkgo_seqstat,
  Quercus = quercus_seqstat,
  Salix = salix_seqstat
)

# Flatten the nested list into a single list of data frames
flattened_seqstat_data <- do.call(c, seqstat_data)

# Combine all data frames into a single data frame
combined_seqstat_data <- bind_rows(flattened_seqstat_data, .id = "Source")

# View the combined data
print(head(combined_seqstat_data))  # Preview the first few rows

# read checkm2 ########################################################################################

# Function to fetch and read .txt files into a list of data frames
read_tsv_files <- function(directory) {
  # Get all .txt file paths in the directory
  file_paths <- list.files(directory, pattern = "\\.tsv$", full.names = TRUE)
  
  # Read each .txt file into a data frame
  data_list <- lapply(file_paths, function(file) {
    read_tsv(file) 
  })
  
  names(data_list) <- basename(file_paths) # Use file names as list names
  
  return(data_list)
}

# Load data from each directory
betula_checkm <- read_tsv_files(betula_path)
ginkgo_checkm <- read_tsv_files(ginkgo_path)
quercus_checkm <- read_tsv_files(quercus_path)
salix_checkm <- read_tsv_files(salix_path)

# Combine all lists into a single list of data frames
checkm_data <- list(
  Betula = betula_checkm,
  Ginkgo = ginkgo_checkm,
  Quercus = quercus_checkm,
  Salix = salix_checkm
)

# Flatten the nested list into a single list of data frames
flattened_checkm_data <- do.call(c, checkm_data)

# Combine all data frames into a single data frame
combined_checkm_data <- bind_rows(flattened_checkm_data, .id = "Source")

# View the combined data
print(head(combined_checkm_data))  # Preview the first few rows

# prepare the dataframes
combined_seqstat_data <- combined_seqstat_data %>%
  rename("Name" = "Bin:",
         "Source_Seqstat" = "Source")

#join the data
Dat <- combined_checkm_data %>%
  left_join(combined_seqstat_data, by = "Name")

# export data
write.csv(Dat, file = "Raw_Data/Final_Bins.csv", row.names = FALSE)

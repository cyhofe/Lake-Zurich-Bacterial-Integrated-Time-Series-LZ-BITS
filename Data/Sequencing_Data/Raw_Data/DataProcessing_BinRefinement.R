#clear R's brain
rm(list = ls())

#libraries
library(utils)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)

#set working directory
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Sequencing_Data/")

#load data
betula_path <- "Raw_Data/Betula_Download/4_Binning/Bin_Refinement/"
ginkgo_path <- "Raw_Data/Ginkgo_Download/4_Binning/Bin_Refinement/"
quercus_path <- "Raw_Data/Quercus_Download/4_Binning/Bin_Refinement/"
salix_path <- "Raw_Data/Salix_Download/4_Binning/Bin_Refinement/"

# Function to fetch and read .tsv files into a list of data frames
read_tsv_files <- function(directory) {
  # Get all .tsv file paths in the directory
  file_paths <- list.files(directory, pattern = "\\.tsv$", full.names = TRUE)
  
  # Read each .tsv file into a data frame and store in a named list
  data_list <- lapply(file_paths, read_tsv)
  names(data_list) <- basename(file_paths) # Use file names as list names
  
  return(data_list)
}

# Load data from each directory
betula_data <- read_tsv_files(betula_path)
ginkgo_data <- read_tsv_files(ginkgo_path)
quercus_data <- read_tsv_files(quercus_path)
salix_data <- read_tsv_files(salix_path)

# Combine all lists into a single list of data frames
all_data <- list(
  Betula = betula_data,
  Ginkgo = ginkgo_data,
  Quercus = quercus_data,
  Salix = salix_data
)

# Filter for data frames containing '_all_bins_' in their names
all_bins_data <- lapply(all_data, function(data_list) {
  data_list[grep("_all_bins_", names(data_list))]
})

# Combine filtered data frames into a single list
combined_all_bins_data <- unlist(all_bins_data, recursive = FALSE)

# Filter for data frames containing '_all_bins_' in their names
final_bins_data <- lapply(all_data, function(data_list) {
  data_list[grep("_final_bins_", names(data_list))]
})

# Combine filtered data frames into a single list
combined_final_bins_data <- unlist(final_bins_data, recursive = FALSE)

# Combine all data frames from the list into one data frame
combined_all_bins_data <- bind_rows(combined_all_bins_data, .id = "Source")
combined_final_bins_data <- bind_rows(combined_final_bins_data, .id = "Source")

#export data
write.csv(combined_all_bins_data, file = "Raw_Data/BinRefinement_AllBins_Report.csv", row.names = FALSE)
write.csv(combined_final_bins_data, file = "Raw_Data/BinRefinement_RefinedBins_Report.csv", row.names = FALSE)


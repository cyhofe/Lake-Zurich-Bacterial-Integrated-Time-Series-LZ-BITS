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
betula_path <- "Raw_Data/Betula_Download/4_Binning/Bin_Curation/"
ginkgo_path <- "Raw_Data/Ginkgo_Download/4_Binning/Bin_Curation/"
quercus_path <- "Raw_Data/Quercus_Download/4_Binning/Bin_Curation/"
salix_path <- "Raw_Data/Salix_Download/4_Binning/Bin_Curation/"

# Function to fetch and read .tsv files into a list of data frames
read_tsv_files <- function(directory) {
  # Get all .tsv file paths in the directory
  file_paths <- list.files(directory, pattern = "\\.tsv$", full.names = TRUE)
  
  # Read each .tsv file into a data frame and store in a named list
  data_list <- lapply(file_paths, function(file) {
    read_tsv(file, col_names = c("bin_id", "completeness", "contamination")) # Explicit column names
  })
  
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

# High-quality bins
highquality_bins_data <- lapply(all_data, function(data_list) {
  data_list[grep("HighQuality", names(data_list))]
})

# Combine HighQuality bins into one data frame
combined_highquality_bins_data <- bind_rows(
  unlist(highquality_bins_data, recursive = FALSE),
  .id = "Source"
)

# Low-quality bins
lowquality_bins_data <- lapply(all_data, function(data_list) {
  data_list[grep("LowQuality", names(data_list))]
})

# Combine LowQuality bins into one data frame
combined_lowquality_bins_data <- bind_rows(
  unlist(lowquality_bins_data, recursive = FALSE),
  .id = "Source"
)

# Removed bins
removed_bins_data <- lapply(all_data, function(data_list) {
  data_list[grep("Removed", names(data_list))]
})

# Combine Removed bins into one data frame
combined_removed_bins_data <- bind_rows(
  unlist(removed_bins_data, recursive = FALSE),
  .id = "Source"
)

# View combined data frames
print(combined_highquality_bins_data)
print(combined_lowquality_bins_data)
print(combined_removed_bins_data)

# combine data
Dat <- rbind(combined_highquality_bins_data, combined_lowquality_bins_data, combined_removed_bins_data)

#add missing variables
Dat <- Dat %>%
  mutate(BinType = str_extract(Source, "HighQuality|LowQuality|Removed"),
         Sample = str_extract(Source, "NZE\\d{2}[A-Z]{3}\\d{4}"),
         Date = str_extract(Source, "\\d{2}[A-Z]{3}\\d{2}"),
         Date = dmy(paste0(Date)))

#export the data
write.csv(Dat, file = "Raw_Data/BinCuration_Report.csv", row.names = FALSE)
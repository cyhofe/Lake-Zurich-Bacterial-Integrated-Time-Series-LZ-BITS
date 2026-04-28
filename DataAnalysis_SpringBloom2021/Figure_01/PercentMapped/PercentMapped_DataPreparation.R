################################################################################
#                                                                              #
#                   Percent Mapped                                             #
#                                                                              #
#   Author  : Cyrill Hofer                                                     #
#   Date    : 09.04.2025                                                       #
#   Purpose : Data Preparation                                                 #
#                                                                              #
################################################################################

################################################################################
#                              CLEAR R'S ENVIRONMENT                           #
################################################################################
# Clear all objects from R's workspace
rm(list = ls())


################################################################################
#                              LOAD LIBRARIES                                  #
################################################################################
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)


################################################################################
#                              SET WORKING DIRECTORY                           #
################################################################################
# Adjust the path as per your project location
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Method_Paper_Analysis/FinalAnalysis/Figure_01/PercentMapped/")


################################################################################
#                              LOAD DATASETS                                   #
################################################################################
mapped_dat <- read_delim("../../../../Data/Species_Data/Raw_Data/MergedMappingStatistics.txt", col_names = F)

################################################################################
#                              PEPARE DATASETS                                 #
################################################################################
mapped_dat <- mapped_dat %>%
  # Extract the date string (e.g., "08MAR21") from the first column
  mutate(Date = str_extract(X1, "\\d{2}[A-Z]{3}\\d{2}"),
         # Convert to Date object, then format to "dd-mm-yyyy"
         Date = dmy(Date)) %>%
  # Rename columns
  rename(
    Mapped_Reads = X2,
    Mapped_Percent = X3
  ) %>%
  # Reorder columns
  select(Date, Mapped_Reads, Mapped_Percent)

################################################################################
#                              EXPORT DATASET                                  #
################################################################################
# Write the data to csv
write_csv(mapped_dat, "PercentMapped.csv") # Export combined dataset to CSV

#clear R's brain
rm(list = ls())

#libraries
library(utils)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)

# set working directory
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Sequencing_Data/Final_Data/")

# load the data 
# Path to folder
path <- "../Raw_Data/Assembly_Comparison_Statistics/AssemblyStats/"

# Get all .tsv files in folder
files <- list.files(path, pattern = "\\.tsv$", full.names = TRUE)

# Read each .tsv into a dataframe and store in a named list
df_list <- lapply(files, read.delim)

# Assign names to list elements based on file names (without extension)
names(df_list) <- tools::file_path_sans_ext(basename(files))

# See names of all dataframes
names(df_list)

# Access a specific dataframe by file name
df_list[["Nanopore_Raw_Assembly_NZE26APR2102.QC.qc_summary"]]

# Combine all dataframes into one
df_all <- bind_rows(df_list, .id = "Source")

# clean the dataframe
df_clean <- df_all %>%
  mutate(
    # Extract Sample (e.g. "IZE01APR2102")
    Sample = str_extract(Source, "[A-Z]{2}[0-9]{2}[A-Z]{3}[0-9]{4}"),
    
    # Extract AssemblyType
    AssemblyType = case_when(
      str_detect(Source, "Illumina_Raw") ~ "Illumina_Raw",
      str_detect(Source, "Illumina_Filtered") ~ "Illumina_Filtered",
      str_detect(Source, "Nanopore_Raw") ~ "Nanopore_Raw",
      str_detect(Source, "Nanopore_Filtered") ~ "Nanopore_Filtered",
      str_detect(Source, "Nanopore_Corrected") ~ "Nanopore_Corrected",
      TRUE ~ "Other"
    )
  ) %>%
  # Drop Source and prefix since no longer needed
  select(-Source, -prefix)

# read the metaquast input
path2 <- "../Raw_Data/Assembly_Comparison_Statistics/QuastStats/"

# Get all .tsv files in that folder
files2 <- list.files(path2, pattern = "\\.tsv$", full.names = TRUE)

# Read each .tsv into a dataframe and store in a named list
df_list2 <- lapply(files2, read.delim, check.names = FALSE)

# Name the list elements after the filenames (no extension)
names(df_list2) <- tools::file_path_sans_ext(basename(files2))

# Combine all dataframes into one
df_all2 <- bind_rows(df_list2, .id = "Source")

# clean the file
df_all2_clean <- df_all2 %>%
  mutate(
    # take everything before the first dot as Sample
    Sample = str_remove(Source, "\\..*$")
  ) %>%
  rename(AssemblyType = Assembly) %>%
  select(Sample, AssemblyType, everything(), -Source) %>%
  select(N50,L50,`Total length`,`# contigs`, Sample, AssemblyType) %>%
  rename(n_Contigs = `# contigs`, Total_Length = `Total length`)

#combine data
dat <- df_clean %>%
  left_join(df_all2_clean, by = c("AssemblyType", "Sample"))

# export csv
write_csv(dat, file = "AssemblyStats.csv")
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
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Species_Data/")

#load data
Species_Dat <- read.delim(file = "Raw_Data/RepresentativeMAGs.list.txt", header = FALSE, col.names = "MAG_ID")
Checkm2_Dat <- read_tsv(file = "Raw_Data/RepresentativeMagsQualityTable.txt")
GenomicFeatures_Dat <- read.delim(file = "Raw_Data/RepresentativeMagsFeatureTable.txt")
Taxonomy_Dat <- read_tsv(file = "Raw_Data/RepresentativeMAGs_Full_TaxonomyTable.tsv")
rRNA_Dat <- read_delim(file = "Raw_Data/Representative_MAGs.rRNA-SeqReport.txt", col_names = FALSE)
Seqstat_Dat <- read_delim(file = "Raw_Data/RepresentativeMagsSizeTable.Final.txt")
Abundance_Dat <- read_delim(file = "Raw_Data/MergedAbundanceFile.txt")

# load metaquast data
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Sequencing_Data/Final_Data/")

# Path to your folder
folder <- "../Raw_Data/Metaqust_Assembly_Comparison"

# List all report files
files <- list.files(folder, pattern = "\\.metaquast\\.report\\.tsv$", full.names = TRUE)

# Function to read one file
read_metaquast <- function(file) {
  # Extract sample name (remove path + extension)
  sample_name <- basename(file) %>%
    str_remove("\\.metaquast\\.report\\.tsv$")
  
  # Extract date string (characters 3–9, e.g. "01APR21")
  date_str <- substr(sample_name, 3, 9)
  
  # Parse to proper date (01APR21 → 01-04-2021)
  date <- format(as.Date(date_str, format = "%d%b%y"), "%d-%m-%Y")
  
  # Read the file
  df <- read_tsv(file, show_col_types = FALSE)
  
  # Add sample + date columns
  df$Sample <- sample_name
  df$Date <- date
  
  return(df)
}

# Read all into a named list
dfs <- lapply(files, read_metaquast)
names(dfs) <- basename(files) %>%
  str_remove("\\.metaquast\\.report\\.tsv$")

combined <- dfs %>%
  bind_rows() %>%
  rename(
    sample                      = Sample,
    date                        = Date,
    assembly                    = Assembly,
    `total length`              = `Total length`,
    `# contigs`                 = `# contigs`,
    n50                         = N50,
    `mismatches per 100 kbp`    = `# mismatches per 100 kbp`,
    `indels per 100 kbp`        = `# indels per 100 kbp`,
    `Total aligned length` = `Total aligned length`) %>%
  select(sample, date, assembly, `total length`, `# contigs`, n50,
         `mismatches per 100 kbp`, `indels per 100 kbp`, `Total aligned length`)


write_csv(combined, "metaquast_summary.csv")
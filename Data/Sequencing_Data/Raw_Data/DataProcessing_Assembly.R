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
betula_dat <- read_delim(file = "Raw_Data/Betula_Download/2_Assembly/Assembly_Report.txt")
ginkgo_dat <- read_delim(file = "Raw_Data/Ginkgo_Download/2_Assembly/Assembly_Report.txt")
quercus_dat <- read_delim(file = "Raw_Data/Quercus_Download/2_Assembly/Assembly_Report.txt")
salix_dat <- read_delim(file = "Raw_Data/Salix_Download/2_Assembly/Assembly_Report.txt")

# prepare the data
Dat <- rbind(betula_dat, ginkgo_dat, quercus_dat, salix_dat)

# add and format the date
Dat <- Dat %>%
  mutate(
    extracted_date = str_extract(Sample, "\\d{2}[A-Z]{3}\\d{2}"), # Extract date pattern
    Date = dmy(paste0(extracted_date)) # Add the year prefix and convert to Date
  ) %>%
  select(-extracted_date) %>%
  mutate(
    Data_Type = if_else(Data_Type == "contigs", "FilteredContigs", "AllContigs") # Correct usage of if_else
  )

# export the data
write.csv(Dat, file = "Raw_Data/Assembly_Report.csv", row.names = FALSE)

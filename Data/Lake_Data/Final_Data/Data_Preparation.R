#clear R's brain
rm(list = ls())

#libraries
library(utils)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(readxl)

#set working directory
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Lake_Data/")

#load data
ST_data_abiotic <- read_excel("Raw_Data/Data_2021.xlsx", sheet = "abiotic")
ST_data_biotic <- read_excel("Raw_Data/Data_2021.xlsx", sheet = "biotic")
LT_data_abiotic <- read_excel("Raw_Data/Data_2009-2023.xlsx", sheet = "abiotic")

# change names
ST_data_abiotic <- ST_data_abiotic %>%
  rename("Date" = "date",
         "Depth" = "depth",
         "Temperature" = "temperature",
         "Conductivity" = "conductivity",
         "Oxygen" = "oxygen",
         "Turbidity" = "turbidity",
         "Chlorophyll" = "chlorophyll")

ST_data_biotic <- ST_data_biotic %>%
  rename("Date" = "date",
         "Depth" = "depth",
         "HNF_CellsPerMl_Epilimnion" = "hnf",
         "Bacteria_CellsPerMl_Epilimnion" = "bacteria")

LT_data_abiotic <- LT_data_abiotic %>%
  rename("Date" = "date",
         "Depth" = "depth",
         "Oxygen" = "oxygen",
         "Phosphate" = "phosphate",
         "Temperature" = "temperature")

# combine short and long term
SpringBloom2021_Metadata <- ST_data_abiotic %>%
  left_join(ST_data_biotic %>% select(-Depth), by = "Date")

# export the data
write.csv(SpringBloom2021_Metadata, file = "Final_Data/SpringBloom2021_Metadata.csv", row.names = FALSE)

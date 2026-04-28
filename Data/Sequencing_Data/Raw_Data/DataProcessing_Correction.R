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

################## Medaka 1 ############################################################################################################

#load data
betula_dat <- read_delim(file = "Raw_Data/Betula_Download/3_Correction/1_Medaka/Correction_1_Medaka_QualityReport.txt")
ginkgo_dat <- read_delim(file = "Raw_Data/Ginkgo_Download/3_Correction/1_Medaka/Correction_1_Medaka_QualityReport.txt")
quercus_dat <- read_delim(file = "Raw_Data/Quercus_Download/3_Correction/1_Medaka/Correction_1_Medaka_QualityReport.txt")
salix_dat <- read_delim(file = "Raw_Data/Salix_Download/3_Correction/1_Medaka/Correction_1_Medaka_QualityReport.txt")

# prepare the data
Dat1 <- rbind(betula_dat, ginkgo_dat, quercus_dat, salix_dat)

# add and format the date
Dat1 <- Dat1 %>%
  mutate(
    extracted_date = str_extract(Sample, "\\d{2}[A-Z]{3}\\d{2}"), # Extract date pattern
    Date = dmy(paste0(extracted_date)) # Add the year prefix and convert to Date
  ) %>%
  select(-extracted_date) %>%
  mutate(Correction_Step = "1",
         Correction_Method = "Medaka",
         Total_Changes = NA)

# change total mapped reads from 0 to na -> missing values, not that 0 mapped
Dat1$Mapped_Reads[Dat1$Mapped_Reads == 0] <- NA

################## Medaka 2 ############################################################################################################

#load data
betula_dat <- read_delim(file = "Raw_Data/Betula_Download/3_Correction/2_Medaka/Correction_2_Medaka_QualityReport.txt")
ginkgo_dat <- read_delim(file = "Raw_Data/Ginkgo_Download/3_Correction/2_Medaka/Correction_2_Medaka_QualityReport.txt")
quercus_dat <- read_delim(file = "Raw_Data/Quercus_Download/3_Correction/2_Medaka/Correction_2_Medaka_QualityReport.txt")
salix_dat <- read_delim(file = "Raw_Data/Salix_Download/3_Correction/2_Medaka/Correction_2_Medaka_QualityReport.txt")

# prepare the data
Dat2 <- rbind(betula_dat, ginkgo_dat, quercus_dat, salix_dat)

# add and format the date
Dat2 <- Dat2 %>%
  mutate(
    extracted_date = str_extract(Sample, "\\d{2}[A-Z]{3}\\d{2}"), # Extract date pattern
    Date = dmy(paste0(extracted_date)) # Add the year prefix and convert to Date
  ) %>%
  select(-extracted_date) %>%
  mutate(Correction_Step = "2",
         Correction_Method = "Medaka",
         Total_Changes = NA)

################## Pilon 1 ############################################################################################################

#load data
betula_dat <- read_delim(file = "Raw_Data/Betula_Download/3_Correction/3_Pilon/Correction_3_Pilon_QualityReport.txt")
ginkgo_dat <- read_delim(file = "Raw_Data/Ginkgo_Download/3_Correction/3_Pilon/Correction_3_Pilon_QualityReport.txt")
quercus_dat <- read_delim(file = "Raw_Data/Quercus_Download/3_Correction/3_Pilon/Correction_3_Pilon_QualityReport.txt")
salix_dat <- read_delim(file = "Raw_Data/Salix_Download/3_Correction/3_Pilon/Correction_3_Pilon_QualityReport.txt")

# prepare the data
Dat3 <- rbind(betula_dat, ginkgo_dat, quercus_dat, salix_dat)

# add and format the date
Dat3 <- Dat3 %>%
  mutate(
    extracted_date = str_extract(Sample, "\\d{2}[A-Z]{3}\\d{2}"), # Extract date pattern
    Date = dmy(paste0(extracted_date)) # Add the year prefix and convert to Date
  ) %>%
  select(-extracted_date) %>%
  mutate(Correction_Step = "3",
         Correction_Method = "Pilon")

################## Pilon 2 ############################################################################################################

#load data
betula_dat <- read_delim(file = "Raw_Data/Betula_Download/3_Correction/4_Pilon/Correction_4_Pilon_QualityReport.txt")
ginkgo_dat <- read_delim(file = "Raw_Data/Ginkgo_Download/3_Correction/4_Pilon/Correction_4_Pilon_QualityReport.txt")
quercus_dat <- read_delim(file = "Raw_Data/Quercus_Download/3_Correction/4_Pilon/Correction_4_Pilon_QualityReport.txt")
salix_dat <- read_delim(file = "Raw_Data/Salix_Download/3_Correction/4_Pilon/Correction_4_Pilon_QualityReport.txt")

# prepare the data
Dat4 <- rbind(betula_dat, ginkgo_dat, quercus_dat, salix_dat)

# add and format the date
Dat4 <- Dat4 %>%
  mutate(
    extracted_date = str_extract(Sample, "\\d{2}[A-Z]{3}\\d{2}"), # Extract date pattern
    Date = dmy(paste0(extracted_date)) # Add the year prefix and convert to Date
  ) %>%
  select(-extracted_date) %>%
  mutate(Correction_Step = "4",
         Correction_Method = "Pilon")

################## Pilon 3 ############################################################################################################

#load data
betula_dat <- read_delim(file = "Raw_Data/Betula_Download/3_Correction/5_Pilon/Correction_5_Pilon_QualityReport.txt")
ginkgo_dat <- read_delim(file = "Raw_Data/Ginkgo_Download/3_Correction/5_Pilon/Correction_5_Pilon_QualityReport.txt")
quercus_dat <- read_delim(file = "Raw_Data/Quercus_Download/3_Correction/5_Pilon/Correction_5_Pilon_QualityReport.txt")
salix_dat <- read_delim(file = "Raw_Data/Salix_Download/3_Correction/5_Pilon/Correction_5_Pilon_QualityReport.txt")

# prepare the data
Dat5 <- rbind(betula_dat, ginkgo_dat, quercus_dat, salix_dat)

# add and format the date
Dat5 <- Dat5 %>%
  mutate(
    extracted_date = str_extract(Sample, "\\d{2}[A-Z]{3}\\d{2}"), # Extract date pattern
    Date = dmy(paste0(extracted_date)) # Add the year prefix and convert to Date
  ) %>%
  select(-extracted_date) %>%
  mutate(Correction_Step = "5",
         Correction_Method = "Pilon")

# combine the data
Dat <- rbind(Dat1, Dat2, Dat3, Dat4, Dat5)

# export the data
write.csv(Dat, file = "Raw_Data/Correction_Report.csv", row.names = FALSE)

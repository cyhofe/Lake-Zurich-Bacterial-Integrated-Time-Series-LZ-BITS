#clear R's brain
rm(list = ls())

#libraries
library(utils)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(readxl)
library(lubridate)
library(data.table)

#set working directory
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Proteomics_Data/Final_Data/74Samples_uniprotkb_taxonomy_id_2_AND_reviewed_tr_2024_12_20_ContaminantsAdded_Database_Run11FEB25/")

#load data
Peptide_Table <- read_excel(path = "/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Proteomics_Data/Raw_Data/74Samples_uniprotkb_taxonomy_id_2_AND_reviewed_tr_2024_12_20_ContaminantsAdded_Database_Run11FEB25/proteinAbundances_Mock_Prepared.xlsx", sheet = 1)
Protein_Table_MedPolish <- read_excel(path = "/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Proteomics_Data/Raw_Data/74Samples_uniprotkb_taxonomy_id_2_AND_reviewed_tr_2024_12_20_ContaminantsAdded_Database_Run11FEB25/proteinAbundances_Mock_Prepared.xlsx", sheet = 3)
Protein_Table_IBAQ <- read_excel(path = "/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Proteomics_Data/Raw_Data/74Samples_uniprotkb_taxonomy_id_2_AND_reviewed_tr_2024_12_20_ContaminantsAdded_Database_Run11FEB25/proteinAbundances_Mock_Prepared.xlsx", sheet = 4)

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
# Peptide table

# Select and rename columns
Peptide_Table <- Peptide_Table %>%
  select(protein_Id, peptide_Id, IDcolumn, nrPeptides, c(10:84)) %>%
  rename("Protein_Group" = "protein_Id",
         "Peptide_Sequence" = "peptide_Id",
         "Protein_ID" = "IDcolumn",
         "Nr_Peptides" = "nrPeptides")

# Add missing column MAG_ID and Peptide_ID
Peptide_Table <- Peptide_Table %>%
  mutate(Species_ID = str_split_fixed(Protein_ID, "-", 2)[, 1]) %>%
  group_by(Protein_Group) %>%
  mutate(Peptide_ID = paste0(Protein_ID, "_p", row_number())) %>%
  ungroup()

# Reorder the columns
Peptide_Table <- Peptide_Table %>%
  select(Species_ID, Protein_ID, Protein_Group, Peptide_ID, Peptide_Sequence, Nr_Peptides, everything())

# Transform to long format
Peptide_Table <- Peptide_Table %>%
  pivot_longer(
    cols = starts_with("ZE-"),       # Select columns that start with "ZE-"
    names_to = "Sample",               # New column for the original column names
    values_to = "Peptide_Abundance"  # New column for the abundance values
  )

# Add date column
Peptide_Table <- Peptide_Table %>%
  mutate(
    Date = str_extract(Sample, "\\d{2}[A-Z]{3}\\d{2}"),  # Extract the date pattern
    Date = dmy(Date)  # Convert the extracted string into Date format
  )

# Generate a Presence column
Presence_Data_Peptide_Table <- Peptide_Table %>%
  group_by(Peptide_ID, Date) %>%
  summarize(
    Peptide_Presence = c("none", "one", "both")[min(sum(!is.na(Peptide_Abundance)), 2) + 1],
    .groups = "drop"
  )

# compute the mean if both present and otherwise take the value present in one
Abundance_Data_Peptide_Table <- Peptide_Table %>%
  group_by(Peptide_ID, Date) %>%
  mutate(Peptide_Abundance = mean(Peptide_Abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-Sample) %>%
  distinct()

# combine the data
Peptide_Table_Combined <- Abundance_Data_Peptide_Table %>%
  left_join(Presence_Data_Peptide_Table, by = c("Peptide_ID", "Date")) %>% 
  mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))

# rearrange columns
Peptide_Table_Combined <- Peptide_Table_Combined %>%
  select(Date, everything()) %>%
  mutate(Date = as.Date(Date))

# calculate absolute and realative amount of peptide per sample
Peptide_Table_Combined <- Peptide_Table_Combined %>%
  group_by(Date) %>%
  mutate(
    Total_Peptide_Abundance_Per_Date = sum(Peptide_Abundance, na.rm = TRUE),
    Relative_Peptide_Abundance_Per_Date = Peptide_Abundance / Total_Peptide_Abundance_Per_Date) %>%
  ungroup()

# filter peptides that solely identify a protein
Peptide_Table_Combined_Filtered <- Peptide_Table_Combined %>%
  filter(Nr_Peptides >= 2) %>%
  distinct() %>%
  mutate(Date = as.Date(Date))

# sanity check
Summary_Peptide_Table <- Peptide_Table %>%
  summarize(Unique_Peptide_IDs = n_distinct(Peptide_ID),
            Unique_Protein_IDs = n_distinct(Protein_ID))

Summary_Peptide_Table_Combined <- Peptide_Table_Combined %>%
  summarize(Unique_Peptide_IDs = n_distinct(Peptide_ID),
            Unique_Protein_IDs = n_distinct(Protein_ID))

Summary_Peptide_Table_Combined_Filtered <- Peptide_Table_Combined_Filtered %>%
  summarize(Unique_Peptide_IDs = n_distinct(Peptide_ID),
            Unique_Protein_IDs = n_distinct(Protein_ID))

# export the datasets
write_csv(Peptide_Table,"Peptide_Table_Replicates.csv")
write_csv(Peptide_Table_Combined, "Peptide_Table_Combined.csv")
write_csv(Peptide_Table_Combined_Filtered,"Peptide_Table_Combined_Filtered.csv")

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
# Protein Table IBAQ

# Select and rename columns
Protein_Table_IBAQ <- Protein_Table_IBAQ %>%
  select(-description, -isotopeLabel, -fasta.id, -c(81:152)) %>%
  rename("Protein_Group" = "protein_Id",
         "Protein_Length" = "protein_length",
         "Protein_ID" = "IDcolumn",
         "Nr_Peptides" = "nrPeptides",
         "Nr_Tryptic_Peptides" = "nr_tryptic_peptides")

# Add missing column MAG_ID
Protein_Table_IBAQ <- Protein_Table_IBAQ %>%
  mutate(Species_ID = str_split_fixed(Protein_ID, "-", 2)[, 1])

# Reorder the columns
Protein_Table_IBAQ <- Protein_Table_IBAQ %>%
  select(Species_ID, Protein_ID, Protein_Group, Protein_Length, Nr_Peptides, Nr_Tryptic_Peptides, everything())

# Transform to long format
Protein_Table_IBAQ <- Protein_Table_IBAQ %>%
  pivot_longer(
    cols = starts_with("ZE-"),       # Select columns that start with "ZE-"
    names_to = "Sample",               # New column for the original column names
    values_to = "Protein_Abundance_IBAQ"  # New column for the abundance values
  )  %>%
  mutate(Sample = str_split_fixed(Sample, "_", 2)[, 1])

# Add date column
Protein_Table_IBAQ <- Protein_Table_IBAQ %>%
  mutate(
    Date = str_extract(Sample, "\\d{2}[A-Z]{3}\\d{2}"),  # Extract the date pattern
    Date = dmy(Date)  # Convert the extracted string into Date format
  )

# Add group size column
Protein_Table_IBAQ <- Protein_Table_IBAQ %>%
  mutate(Protein_Group_Size = str_count(Protein_Group, ";") + 1)

# Generate a Presence column
Presence_Data_Protein_Table_IBAQ <- Protein_Table_IBAQ %>%
  group_by(Protein_ID, Date) %>%
  summarize(
    Protein_Presence = c("none", "one", "both")[min(sum(!is.na(Protein_Abundance_IBAQ)), 2) + 1],
    .groups = "drop"
  )

# compute the mean if both present and otherwise take the value present in one
Abundance_Data_Protein_Table_IBAQ <- Protein_Table_IBAQ %>%
  group_by(Protein_ID, Date) %>%
  mutate(Protein_Abundance_IBAQ = mean(Protein_Abundance_IBAQ, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-Sample) %>%
  distinct()

# combine the data
Protein_Table_IBAQ_Combined <- Abundance_Data_Protein_Table_IBAQ %>%
  left_join(Presence_Data_Protein_Table_IBAQ, by = c("Protein_ID", "Date")) %>% 
  mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))

# calculate absolute and realative amount of peptide per sample
Protein_Table_IBAQ_Combined <- Protein_Table_IBAQ_Combined %>%
  group_by(Date) %>%
  mutate(
    Total_Protein_Abundance_IBAQ_Per_Date = sum(Protein_Abundance_IBAQ, na.rm = TRUE),
    Relative_Protein_Abundance_IBAQ_Per_Date = Protein_Abundance_IBAQ / Total_Protein_Abundance_IBAQ_Per_Date) %>%
  ungroup()

# rearrange columns
Protein_Table_IBAQ_Combined <- Protein_Table_IBAQ_Combined %>%
  select(Date, Species_ID, Protein_ID, Protein_Length, Protein_Group, Protein_Group_Size, Protein_Abundance_IBAQ, Relative_Protein_Abundance_IBAQ_Per_Date, Total_Protein_Abundance_IBAQ_Per_Date, everything())

# filter peptides that solely identify a protein
Protein_Table_IBAQ_Combined_Filtered <- Protein_Table_IBAQ_Combined %>%
  filter(Nr_Peptides >= 2) %>%
  distinct()

# sanity check
Summary_Protein_Table_IBAQ <- Protein_Table_IBAQ %>%
  select(Nr_Peptides, Protein_ID) %>%
  group_by(Protein_ID) %>%
  distinct() %>%
  ungroup() %>%
  summarize(Unique_Peptide_IDs = sum(Nr_Peptides),
            Unique_Protein_IDs = n_distinct(Protein_ID))

Summary_Protein_Table_IBAQ_Combined <- Protein_Table_IBAQ_Combined %>%
  select(Nr_Peptides, Protein_ID) %>%
  group_by(Protein_ID) %>%
  distinct() %>%
  ungroup() %>%
  summarize(Unique_Peptide_IDs = sum(Nr_Peptides),
            Unique_Protein_IDs = n_distinct(Protein_ID))

Summary_Protein_Table_IBAQ_Combined_Filtered <- Protein_Table_IBAQ_Combined_Filtered %>%
  select(Nr_Peptides, Protein_ID) %>%
  group_by(Protein_ID) %>%
  distinct() %>%
  ungroup() %>%
  summarize(Unique_Peptide_IDs = sum(Nr_Peptides),
            Unique_Protein_IDs = n_distinct(Protein_ID))

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
# Protein Table Medpolish

# Select and rename columns
Protein_Table_MedPolish <- Protein_Table_MedPolish %>%
  select(-description, -isotopeLabel, -fasta.id, -c(81:152)) %>%
  rename("Protein_Group" = "protein_Id",
         "Protein_Length" = "protein_length",
         "Protein_ID" = "IDcolumn",
         "Nr_Peptides" = "nrPeptides",
         "Nr_Tryptic_Peptides" = "nr_tryptic_peptides")

# Add missing column MAG_ID
Protein_Table_MedPolish <- Protein_Table_MedPolish %>%
  mutate(Species_ID = str_split_fixed(Protein_ID, "-", 2)[, 1])

# Reorder the columns
Protein_Table_MedPolish <- Protein_Table_MedPolish %>%
  select(Species_ID, Protein_ID, Protein_Group, Protein_Length, Nr_Peptides, Nr_Tryptic_Peptides, everything())

# Transform to long format
Protein_Table_MedPolish <- Protein_Table_MedPolish %>%
  pivot_longer(
    cols = starts_with("ZE-"),       # Select columns that start with "ZE-"
    names_to = "Sample",               # New column for the original column names
    values_to = "Protein_Abundance_MedPolish"  # New column for the abundance values
  ) %>%
  mutate(Sample = str_split_fixed(Sample, "_", 2)[, 1])

# Add date column
Protein_Table_MedPolish <- Protein_Table_MedPolish %>%
  mutate(
    Date = str_extract(Sample, "\\d{2}[A-Z]{3}\\d{2}"),  # Extract the date pattern
    Date = dmy(Date)  # Convert the extracted string into Date format
  )

# Add group size column
Protein_Table_MedPolish <- Protein_Table_MedPolish %>%
  mutate(Protein_Group_Size = str_count(Protein_Group, ";") + 1)

# Generate a Presence column
Presence_Data_Protein_Table_MedPolish <- Protein_Table_MedPolish %>%
  group_by(Protein_ID, Date) %>%
  summarize(
    Protein_Presence = c("none", "one", "both")[min(sum(!is.na(Protein_Abundance_MedPolish)), 2) + 1],
    .groups = "drop"
  )

# compute the mean if both present and otherwise take the value present in one
Abundance_Data_Protein_Table_MedPolish <- Protein_Table_MedPolish %>%
  group_by(Protein_ID, Date) %>%
  mutate(Protein_Abundance_MedPolish = mean(Protein_Abundance_MedPolish, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-Sample) %>%
  distinct()

# combine the data
Protein_Table_MedPolish_Combined <- Abundance_Data_Protein_Table_MedPolish %>%
  left_join(Presence_Data_Protein_Table_MedPolish, by = c("Protein_ID", "Date")) %>% 
  mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))

# calculate absolute and realative amount of peptide per sample
Protein_Table_MedPolish_Combined <- Protein_Table_MedPolish_Combined %>%
  group_by(Date) %>%
  mutate(
    Total_Protein_Abundance_MedPolish_Per_Date = sum(Protein_Abundance_MedPolish, na.rm = TRUE),
    Relative_Protein_Abundance_MedPolish_Per_Date = Protein_Abundance_MedPolish / Total_Protein_Abundance_MedPolish_Per_Date) %>%
  ungroup()

# rearrange columns
Protein_Table_MedPolish_Combined <- Protein_Table_MedPolish_Combined %>%
  select(Date, Species_ID, Protein_ID, Protein_Length, Protein_Group, Protein_Group_Size, Protein_Abundance_MedPolish, Relative_Protein_Abundance_MedPolish_Per_Date, Total_Protein_Abundance_MedPolish_Per_Date, everything())

# filter peptides that solely identify a protein
Protein_Table_MedPolish_Combined_Filtered <- Protein_Table_MedPolish_Combined %>%
  filter(Nr_Peptides >= 2) %>%
  distinct()

# sanity check
Summary_Protein_Table_MedPolish <- Protein_Table_MedPolish %>%
  select(Nr_Peptides, Protein_ID) %>%
  group_by(Protein_ID) %>%
  distinct() %>%
  ungroup() %>%
  summarize(Unique_Peptide_IDs = sum(Nr_Peptides),
            Unique_Protein_IDs = n_distinct(Protein_ID))

Summary_Protein_Table_MedPolish_Combined <- Protein_Table_MedPolish_Combined %>%
  select(Nr_Peptides, Protein_ID) %>%
  group_by(Protein_ID) %>%
  distinct() %>%
  ungroup() %>%
  summarize(Unique_Peptide_IDs = sum(Nr_Peptides),
            Unique_Protein_IDs = n_distinct(Protein_ID))

Summary_Protein_Table_MedPolish_Combined_Filtered <- Protein_Table_MedPolish_Combined_Filtered %>%
  select(Nr_Peptides, Protein_ID) %>%
  group_by(Protein_ID) %>%
  distinct() %>%
  ungroup() %>%
  summarize(Unique_Peptide_IDs = sum(Nr_Peptides),
            Unique_Protein_IDs = n_distinct(Protein_ID))

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
# Final Protein Table

Protein_Table <- Protein_Table_IBAQ %>% 
  left_join(Protein_Table_MedPolish %>% select(Protein_ID, Date, Sample, Protein_Abundance_MedPolish), by = c("Protein_ID", "Date", "Sample")) %>%
  select(Date, Species_ID, Protein_ID, Protein_Group, Protein_Group_Size, Protein_Abundance_IBAQ, Protein_Abundance_MedPolish, Protein_Length, Nr_Peptides, Nr_Tryptic_Peptides, Sample, everything())
  
Protein_Table_Combined <- Protein_Table_IBAQ_Combined %>%
  left_join(Protein_Table_MedPolish_Combined %>% select(Protein_ID, Date, Protein_Abundance_MedPolish, Relative_Protein_Abundance_MedPolish_Per_Date, Total_Protein_Abundance_MedPolish_Per_Date), by = c("Protein_ID", "Date")) %>%
  mutate(Date = as.Date(Date)) %>%
  select(Date, Species_ID, Protein_ID, Protein_Group, Protein_Group_Size, Protein_Abundance_IBAQ, Relative_Protein_Abundance_IBAQ_Per_Date, Total_Protein_Abundance_IBAQ_Per_Date, Protein_Abundance_MedPolish, Relative_Protein_Abundance_MedPolish_Per_Date, Total_Protein_Abundance_MedPolish_Per_Date, Protein_Length, Nr_Peptides, Nr_Tryptic_Peptides, everything())

Protein_Table_Combined_Filtered <- Protein_Table_IBAQ_Combined_Filtered %>%
  left_join(Protein_Table_MedPolish_Combined_Filtered %>% select(Protein_ID, Date, Protein_Abundance_MedPolish, Relative_Protein_Abundance_MedPolish_Per_Date, Total_Protein_Abundance_MedPolish_Per_Date), by = c("Protein_ID", "Date")) %>%
  mutate(Date = as.Date(Date)) %>%
  select(Date, Species_ID, Protein_ID, Protein_Group, Protein_Group_Size, Protein_Abundance_IBAQ, Relative_Protein_Abundance_IBAQ_Per_Date, Total_Protein_Abundance_IBAQ_Per_Date, Protein_Abundance_MedPolish, Relative_Protein_Abundance_MedPolish_Per_Date, Total_Protein_Abundance_MedPolish_Per_Date, Protein_Length, Nr_Peptides, Nr_Tryptic_Peptides, everything())

# export the datasets
write_csv(Protein_Table, "Protein_Table_Replicates.csv")
write_csv(Protein_Table_Combined, "Protein_Table_Combined.csv")
write_csv(Protein_Table_Combined_Filtered, "Protein_Table_Combined_Filtered.csv")

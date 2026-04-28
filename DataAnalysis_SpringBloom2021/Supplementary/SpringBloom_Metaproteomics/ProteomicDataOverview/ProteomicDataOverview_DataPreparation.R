################################################################################
#                                                                              #
#                             Proteome Overview                                #
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


################################################################################
#                              SET WORKING DIRECTORY                           #
################################################################################
# Adjust the path as per your project location
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Supplementary/SpringBloom_Metaproteomics/ProteomicDataOverview/")


################################################################################
#                              LOAD DATASETS                                   #
################################################################################
# Load only the required columns for peptides
Peptide_Dat <- read_csv("../../../../Datasets/72Samples_CostumDB/Integrated_Peptide_Table.csv", col_select = c(Date, Peptide_ID, Species_ID, Protein_ID ,Peptide_Abundance, Peptide_Presence, Phylum))

# Load only the required columns for proteins
Protein_Dat <- read_csv("../../../../Datasets/72Samples_CostumDB/Integrated_Protein_Table.csv", col_select = c(Date, Protein_ID, Species_ID, Nr_Peptides, Nr_Tryptic_Peptides, Protein_Abundance_IBAQ, Protein_Presence, Phylum, Class, Order, Family, Genus, Species, KEGG_ID))

################################################################################
#                              PEPARE DATASETS                                 #
################################################################################
# ----------------------------Peptide Dataset-----------------------------------
# 1. Peptides per Sample per Dataset (CostumDB)
Peptide_Summary <- Peptide_Dat %>% # Start with CostumDB peptides dataset
  filter(Peptide_Abundance > 0) %>% # Keep only rows with abundance > 0
  mutate(Overall_Unique_Peptides_Detected = n_distinct(Peptide_ID)) %>% # calculate total unique peptides
  group_by(Date) %>% # Group by sample Date
  mutate(                                       
    Unique_Peptides = n_distinct(Peptide_ID), # Count unique peptides per sample
    Total_Peptide_Abundance = sum(Peptide_Abundance, na.rm = TRUE), # Sum total peptide abundance per sample
    Unique_Peptides_In_One_Replicate = sum(Peptide_Presence == "one"), # Count peptides found in one replicate
    Unique_Peptides_In_Both_Replicates = sum(Peptide_Presence == "both"), # Count peptides found in both replicates
    Percent_Peptides_In_One_Replicate = (Unique_Peptides_In_One_Replicate / Unique_Peptides) * 100, # Calculate % of peptides in one replicate
    Percent_Peptides_In_Both_Replicates = (Unique_Peptides_In_Both_Replicates / Unique_Peptides) * 100, # Calculate % of peptides in both replicates
    Abundance_Peptides_In_One_Replicate = sum(ifelse(Peptide_Presence == "one", Peptide_Abundance, 0), na.rm = TRUE), # Sum abundance of peptides in one replicate
    Abundance_Peptides_In_Both_Replicates = sum(ifelse(Peptide_Presence == "both", Peptide_Abundance, 0), na.rm = TRUE), # Sum abundance of peptides in both replicates
    Percent_Abundance_Peptides_In_One_Replicate = (Abundance_Peptides_In_One_Replicate / Total_Peptide_Abundance) * 100, # % abundance for one replicate peptides 
    Percent_Abundance_Peptides_In_Both_Replicates = (Abundance_Peptides_In_Both_Replicates / Total_Peptide_Abundance) * 100) %>% # % abundance for both replicate peptides %>%
ungroup() # ungroup the summary per date

# ----------------------------Protein Dataset-----------------------------------
# 1. Proteins per Sample per Dataset (CostumDB)
Protein_Summary <- Protein_Dat %>% # Start with CostumDB proteins dataset
  filter(Protein_Abundance_IBAQ > 0) %>% # Keep only rows with abundance > 0
  mutate(Overall_Unique_Proteins_Detected = n_distinct(Protein_ID)) %>% # calculate total unique proteins
  group_by(Date) %>% # Group by sample Date
  mutate(
    Unique_Proteins = n_distinct(Protein_ID), # Count unique proteins per sample
    Total_Protein_Abundance = sum(Protein_Abundance_IBAQ, na.rm = TRUE), # Sum total protein abundance per sample
    Unique_Proteins_In_One_Replicate = sum(Protein_Presence == "one"), # Count proteins found in one replicate
    Unique_Proteins_In_Both_Replicates = sum(Protein_Presence == "both"), # Count proteins found in both replicates
    Percent_Proteins_In_One_Replicate = (Unique_Proteins_In_One_Replicate / Unique_Proteins) * 100, # Calculate % of proteins in one replicate
    Percent_Proteins_In_Both_Replicates = (Unique_Proteins_In_Both_Replicates / Unique_Proteins) * 100, # Calculate % of proteins in both replicates
    Abundance_Proteins_In_One_Replicate = sum(ifelse(Protein_Presence == "one", Protein_Abundance_IBAQ, 0), na.rm = TRUE), # Sum abundance of proteins in one replicate
    Abundance_Proteins_In_Both_Replicates = sum(ifelse(Protein_Presence == "both", Protein_Abundance_IBAQ, 0), na.rm = TRUE), # Sum abundance of proteins in both replicates
    Percent_Abundance_Proteins_In_One_Replicate = (Abundance_Proteins_In_One_Replicate / Total_Protein_Abundance) * 100, # % abundance for one replicate proteins
    Percent_Abundance_Proteins_In_Both_Replicates = (Abundance_Proteins_In_Both_Replicates / Total_Protein_Abundance) * 100) %>% # % abundance for both replicate proteins
  ungroup() # ungroup the summary by date

################################################################################
#                              EXPORT DATASET                                  #
################################################################################
# Merge the Datasets
Combined_Data <- Peptide_Summary %>%
  left_join(Protein_Summary, by = c("Date", "Species_ID", "Protein_ID", "Phylum"))

# Write the data to csv
write_csv(Combined_Data, "Proteomic_Data_Overview.csv") # Export combined dataset to CSV

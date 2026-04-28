################################################################################
#                                                                              #
#                           Protein Recovered                                  #
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
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Method_Paper_Analysis/FinalAnalysis/Figure_03/ProteinRecovered/")


################################################################################
#                              LOAD DATASETS                                   #
################################################################################
prot_dat <- read_csv("../../../../Datasets/72Samples_CostumDB/Integrated_Protein_Table.csv", col_select = c("Species_ID", "Protein_ID", "Protein_Abundance_IBAQ", "Date", "Species_Abundance", "Phylum"))
gene_dat <- read_csv("../../../../Datasets/72Samples_CostumDB/Integrated_Gene_Table.csv", col_select = c("Species_ID", "Protein_ID"))
################################################################################
#                              PEPARE DATASETS                                 #
################################################################################


# filter species that are not present and fiter contaminants
prot_dat <- prot_dat %>%
  filter(Phylum != "Y-FGCZCont") %>%
  group_by(Species_ID) %>%
  mutate(Total_Species_Abundance = sum(Species_Abundance)) %>%
  filter(Total_Species_Abundance != 0)

# calculate genes per species
genes_per_species <- gene_dat %>%
  group_by(Species_ID) %>%
  summarise(Genes_Per_Species = n_distinct(Protein_ID))

# calculate proteins per species
proteins_per_species <- prot_dat %>%
  filter(Protein_Abundance_IBAQ > 0) %>%
  group_by(Species_ID, Date) %>%
  summarise(Proteins_Per_Species_Per_Date = n_distinct(Protein_ID))

# combine data
combined_dat <- prot_dat %>%
  select(Species_ID, Species_Abundance, Date) %>%
  distinct() %>%
  left_join(proteins_per_species, by = c("Species_ID", "Date")) %>%
  left_join(genes_per_species, by = "Species_ID")

# calculate protein recovery
protein_recovery_summary <- combined_dat %>%
  filter(Proteins_Per_Species_Per_Date > 0) %>%
  group_by(Date) %>%
  summarise(
    Total_Proteins = sum(Proteins_Per_Species_Per_Date),
    Total_Genes = sum(Genes_Per_Species),
    Protein_Recovery_Percent = (Total_Proteins / Total_Genes) * 100
  ) %>%
  ungroup()

################################################################################
#                              EXPORT DATASET                                  #
################################################################################
# Write the data to csv
write_csv(protein_recovery_summary, "ProteinRecovery.csv") # Export combined dataset to CSV

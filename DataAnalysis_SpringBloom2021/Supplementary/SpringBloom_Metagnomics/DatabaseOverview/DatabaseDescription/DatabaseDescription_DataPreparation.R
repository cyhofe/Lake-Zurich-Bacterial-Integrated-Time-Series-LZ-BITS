################################################################################
#                                                                              #
#                             Database Description                             #
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
library(cleaver)
library(purrr)

################################################################################
#                              SET WORKING DIRECTORY                           #
################################################################################
# Adjust the path as per your project location
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Supplementary/SpringBloom_Metagnomics/DatabaseOverview/DatabaseDescription/")


################################################################################
#                              LOAD DATASETS                                   #
################################################################################
gene_dat <- read_csv("../../../../../Data/Gene_Data/Final_Data/Gene_Table.csv")
species_dat <- read_csv("../../../../../Data/Species_Data/Final_Data/Species_Table.csv")

################################################################################
#                              PEPARE DATASETS                                 #
################################################################################
dat <- gene_dat %>%
  rowwise() %>%
  mutate(
    # --- DIA-NN: --met-excision ---
    # Remove the initial methionine (M) if present at the N-terminus
    Cleaned_Protein_Sequence = ifelse(
      str_sub(Protein_Sequence, 1, 1) == "M",
      str_sub(Protein_Sequence, 2),
      Protein_Sequence),
    
    # --- DIA-NN: --cut K*,R* --missed-cleavages 1 ---
    # Perform in silico tryptic digestion using 'cleaver' package
    # Cleave after K or R unless followed by P (handled automatically by cleaver)
    # Allow up to 1 missed cleavage to match DIA-NN setting
    All_Tryptic_Peptides = list(
      cleave(Cleaned_Protein_Sequence, enzym = "trypsin", missedCleavages = 1)),
    
    # --- DIA-NN: --min-pep-len 7 --max-pep-len 30 ---
    # Filter resulting peptides by amino acid length (7 to 30 residues)
    Filtered_Tryptic_Peptides = list(
      All_Tryptic_Peptides %>%
        unlist() %>%
        {.[nchar(.) >= 7 & nchar(.) <= 30]})) %>%
  mutate(
    # Count the number of valid peptides per protein
    Nr_Tryptic_Peptides = length(Filtered_Tryptic_Peptides),
    
    # Concatenate all peptides into a single string (for compact representation)
    Tryptic_Peptides = paste(Filtered_Tryptic_Peptides, collapse = ";"),
    
    # Record the length of each peptide in a semicolon-separated string
    Peptide_Length = paste(nchar(Filtered_Tryptic_Peptides), collapse = ";")) %>%
  ungroup() %>%
  
  # Clean up intermediate columns
  select(-All_Tryptic_Peptides, -Filtered_Tryptic_Peptides, -Cleaned_Protein_Sequence) %>%
  
  # Merge in taxonomic metadata
  left_join(species_dat %>%
              select(Species_ID, Phylum, Class, Order, Family, Genus, Species) %>%
              distinct(),
            by = "Species_ID")

################################################################################
#                              EXPORT DATASET                                  #
################################################################################
# Write the data to csv
write_csv(dat, "ProteinDataAnnotated.csv") # Export combined dataset to CSV


################################################################################
#                                                                              #
#                       Proteome Quality Assessment                            #
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
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Supplementary/SpringBloom_Metaproteomics/ProteinQualityAssessment/")


################################################################################
#                              LOAD DATASETS                                   #
################################################################################
# Load only the required columns for peptides
Peptide_Dat <- read_csv("../../../../Datasets/72Samples_CostumDB/Integrated_Peptide_Table.csv", col_select = c(Date, Peptide_ID, Species_ID, Protein_ID ,Peptide_Abundance))

# Load only the required columns for proteins
Protein_Dat <- read_csv("../../../../Datasets/72Samples_CostumDB/Integrated_Protein_Table.csv", col_select = c(Date, Protein_ID, Species_ID, Protein_Abundance_IBAQ))


################################################################################
#                              PEPARE DATASETS                                 #
################################################################################
dat <- Peptide_Dat %>%
  left_join(Protein_Dat, by = c("Date", "Protein_ID", "Species_ID"))

################################################################################
#                              EXPORT DATASET                                  #
################################################################################
# Write the data to csv
write_csv(dat, "ProteomeQualityAssessment.csv") # Export combined dataset to CSV

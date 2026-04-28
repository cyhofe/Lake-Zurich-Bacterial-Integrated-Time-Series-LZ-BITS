################################################################################
#                                                                              #
#                             Mag Overview                                     #
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


################################################################################
#                              SET WORKING DIRECTORY                           #
################################################################################
# Adjust the path as per your project location
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Method_Paper_Analysis/MagOverview/")


################################################################################
#                              LOAD DATASETS                                   #
################################################################################
dat <- read_csv("../../Data/MAG_Data/Final_Data/MAG_Table.csv")

################################################################################
#                              PEPARE DATASETS                                 #
################################################################################
dat <- dat %>%
select(MAG_ID, Species_ID, Completeness, Contamination, Nr_Contigs, Nr_Bases, Nr_5s, Nr_16s, Nr_23s, Avg_Contig_Size, ANI_To_Reference, Accession_Number, Phylum, Class, Order, Family, Genus, Species, MAGs_Per_Species)

################################################################################
#                              EXPORT DATASET                                  #
################################################################################
# Write the data to csv
write_csv(dat, "MagOverview.csv") # Export combined dataset to CSV

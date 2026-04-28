################################################################################
#                                                                              #
#                             Species O verview                                #
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
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Method_Paper_Analysis/SpeciesOverview/")


################################################################################
#                              LOAD DATASETS                                   #
################################################################################
# load the species dat
dat <- read_csv("../../Data/Species_Data/Final_Data/Species_Table.csv")

################################################################################
#                              PEPARE DATASETS                                 #
################################################################################
dat_all <- dat %>%
  group_by(Species_ID) %>%
  mutate(PerSpecies_Total_Abundance = sum(Species_Abundance)) %>%
  select(Species_ID,Estimated_Genome_Size, GC_Content, Phylum, Class, Order, Family, Genus, Species, ANI_To_Reference, Total_CDS, Coding_Density, PerSpecies_Total_Abundance) %>%
  distinct()

dat_per_date <- dat %>%
  select(Species_ID, Date, Species_Abundance, Phylum)

################################################################################
#                              EXPORT DATASET                                  #
################################################################################
# Write the data to csv
write_csv(dat_all, "SpeciesOverview.csv") # Export combined dataset to CSV
write_csv(dat_per_date, "SpeciesOverview_PerDate.csv")
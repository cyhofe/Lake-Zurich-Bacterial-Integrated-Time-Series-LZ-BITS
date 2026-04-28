################################################################################
#                                                                              #
#                             Database Comparison                              #
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
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Figure_02/DatabaseComparison/")


################################################################################
#                              LOAD DATASETS                                   #
################################################################################
# Load only the required columns for peptides
Peptide_Dat_CostumDB <- read_csv("../../../Datasets/74Samples_CostumDB/Integrated_Peptide_Table.csv",col_select = c(Date, Peptide_ID, Peptide_Abundance, Peptide_Presence, Phylum))
Peptide_Dat_UniProtDB <- read_csv("../../../Datasets/74Samples_UniprotDB/Integrated_Peptide_Table.csv",col_select = c(Date, Peptide_ID, Peptide_Abundance, Peptide_Presence, Phylum))
Peptide_Dat_NcbiRefSeq <- read_csv("../../../Datasets/74Samples_NcbiRefSeq/Integrated_Peptide_Table.csv",col_select = c(Date, Peptide_ID, Peptide_Abundance, Peptide_Presence, Phylum))
Peptide_Dat_MockDB <- read_csv("../../../Datasets/74Samples_Mock/Integrated_Peptide_Table.csv",col_select = c(Date, Peptide_ID, Peptide_Abundance, Peptide_Presence, Phylum))

# Load only the required columns for proteins
Protein_Dat_CostumDB <- read_csv("../../../Datasets/74Samples_CostumDB/Integrated_Protein_Table.csv",col_select = c(Date, Protein_ID, Protein_Abundance_IBAQ, Protein_Presence, Phylum))
Protein_Dat_UniProtDB <- read_csv("../../../Datasets/74Samples_UniprotDB/Integrated_Protein_Table.csv",col_select = c(Date, Protein_ID, Protein_Abundance_IBAQ, Protein_Presence, Phylum))
Protein_Dat_NcbiRefSeq <- read_csv("../../../Datasets/74Samples_NcbiRefSeq/Integrated_Protein_Table.csv",col_select = c(Date, Protein_ID, Protein_Abundance_IBAQ, Protein_Presence, Phylum))
Protein_Dat_MockDB <- read_csv("../../../Datasets/74Samples_Mock/Integrated_Protein_Table.csv",col_select = c(Date, Protein_ID, Protein_Abundance_IBAQ, Protein_Presence, Phylum))

# ---- Peptides ----
Peptide_Dat_CostumDB <- Peptide_Dat_CostumDB %>%
  filter(!grepl("Y-FGCZCont", Peptide_ID))

Peptide_Dat_UniProtDB <- Peptide_Dat_UniProtDB %>%
  filter(!grepl("Cont", Peptide_ID))

Peptide_Dat_NcbiRefSeq <- Peptide_Dat_NcbiRefSeq %>%
  filter(!grepl("Cont", Peptide_ID))

Peptide_Dat_MockDB <- Peptide_Dat_MockDB %>%
  filter(!grepl("Cont", Peptide_ID))


# ---- Proteins ----
Protein_Dat_CostumDB <- Protein_Dat_CostumDB %>%
  filter(!grepl("Y-FGCZCont", Protein_ID))

Protein_Dat_UniProtDB <- Protein_Dat_UniProtDB %>%
  filter(!grepl("Cont", Protein_ID))

Protein_Dat_NcbiRefSeq <- Protein_Dat_NcbiRefSeq %>%
  filter(!grepl("Cont", Protein_ID))

Protein_Dat_MockDB <- Protein_Dat_MockDB %>%
  filter(!grepl("Cont", Protein_ID))


################################################################################
#                              PEPARE DATASETS                                 #
################################################################################
# ----------------------------Peptide Dataset-----------------------------------
# 1. Peptides per Sample per Dataset (CostumDB)
Peptide_Summary_CostumDB <- Peptide_Dat_CostumDB %>% # Start with CostumDB peptides dataset
  filter(Peptide_Abundance > 0) %>% # Keep only rows with abundance > 0
  mutate(Overall_Unique_Peptides_Detected = n_distinct(Peptide_ID)) %>% # count the overall unique peptides
  group_by(Date) %>% # Group by sample Date
  summarize(                                       
    Unique_Peptides = n_distinct(Peptide_ID), # Count unique peptides per sample
    Total_Peptide_Abundance = sum(Peptide_Abundance, na.rm = TRUE), # Sum total peptide abundance per sample
    Unique_Peptides_In_One_Replicate = sum(Peptide_Presence == "one"), # Count peptides found in one replicate
    Unique_Peptides_In_Both_Replicates = sum(Peptide_Presence == "both"), # Count peptides found in both replicates
    Percent_Peptides_In_One_Replicate = (Unique_Peptides_In_One_Replicate / Unique_Peptides) * 100, # Calculate % of peptides in one replicate
    Percent_Peptides_In_Both_Replicates = (Unique_Peptides_In_Both_Replicates / Unique_Peptides) * 100, # Calculate % of peptides in both replicates
    Abundance_Peptides_In_One_Replicate = sum(ifelse(Peptide_Presence == "one", Peptide_Abundance, 0), na.rm = TRUE), # Sum abundance of peptides in one replicate
    Abundance_Peptides_In_Both_Replicates = sum(ifelse(Peptide_Presence == "both", Peptide_Abundance, 0), na.rm = TRUE), # Sum abundance of peptides in both replicates
    Percent_Abundance_Peptides_In_One_Replicate = (Abundance_Peptides_In_One_Replicate / Total_Peptide_Abundance) * 100, # % abundance for one replicate peptides 
    Percent_Abundance_Peptides_In_Both_Replicates = (Abundance_Peptides_In_Both_Replicates / Total_Peptide_Abundance) * 100, # % abundance for both replicate peptides
    Overall_Unique_Peptides_Detected = first(Overall_Unique_Peptides_Detected), # fetch the overall unique peptides detected variable
    .groups = "drop") %>% # Drop grouping after summarizing
  mutate(Database = "CostumDB") # Add a column to indicate the database

# 2. Peptides per Sample per Dataset (UniProtDB)
Peptide_Summary_UniProtDB <- Peptide_Dat_UniProtDB %>% # Start with CostumDB peptides dataset
  filter(Peptide_Abundance > 0) %>% # Keep only rows with abundance > 0
  mutate(Overall_Unique_Peptides_Detected = n_distinct(Peptide_ID)) %>% # calculate the total unique peptides found in the whole dataset
  group_by(Date) %>% # Group by sample Date
  summarize(                                       
    Unique_Peptides = n_distinct(Peptide_ID), # Count unique peptides per sample
    Total_Peptide_Abundance = sum(Peptide_Abundance, na.rm = TRUE), # Sum total peptide abundance per sample
    Unique_Peptides_In_One_Replicate = sum(Peptide_Presence == "one"), # Count peptides found in one replicate
    Unique_Peptides_In_Both_Replicates = sum(Peptide_Presence == "both"), # Count peptides found in both replicates
    Percent_Peptides_In_One_Replicate = (Unique_Peptides_In_One_Replicate / Unique_Peptides) * 100, # Calculate % of peptides in one replicate
    Percent_Peptides_In_Both_Replicates = (Unique_Peptides_In_Both_Replicates / Unique_Peptides) * 100, # Calculate % of peptides in both replicates
    Abundance_Peptides_In_One_Replicate = sum(ifelse(Peptide_Presence == "one", Peptide_Abundance, 0), na.rm = TRUE), # Sum abundance of peptides in one replicate
    Abundance_Peptides_In_Both_Replicates = sum(ifelse(Peptide_Presence == "both", Peptide_Abundance, 0), na.rm = TRUE), # Sum abundance of peptides in both replicates
    Percent_Abundance_Peptides_In_One_Replicate = (Abundance_Peptides_In_One_Replicate / Total_Peptide_Abundance) * 100, # % abundance for one replicate peptides 
    Percent_Abundance_Peptides_In_Both_Replicates = (Abundance_Peptides_In_Both_Replicates / Total_Peptide_Abundance) * 100, # % abundance for both replicate peptides
    Overall_Unique_Peptides_Detected = first(Overall_Unique_Peptides_Detected), # fetch the overall unique peptides detected variable
    .groups = "drop") %>% # Drop grouping after summarizing
  mutate(Database = "UniProtDB") # Add a column to indicate the database

# 3. Peptides per Sample per Dataset (NCBIrefseq)
Peptide_Summary_NcbiRefSeq <- Peptide_Dat_NcbiRefSeq %>% # Start with CostumDB peptides dataset
  filter(Peptide_Abundance > 0) %>% # Keep only rows with abundance > 0
  mutate(Overall_Unique_Peptides_Detected = n_distinct(Peptide_ID)) %>% # calculate the total unique peptides found in the whole dataset
  group_by(Date) %>% # Group by sample Date
  summarize(                                       
    Unique_Peptides = n_distinct(Peptide_ID), # Count unique peptides per sample
    Total_Peptide_Abundance = sum(Peptide_Abundance, na.rm = TRUE), # Sum total peptide abundance per sample
    Unique_Peptides_In_One_Replicate = sum(Peptide_Presence == "one"), # Count peptides found in one replicate
    Unique_Peptides_In_Both_Replicates = sum(Peptide_Presence == "both"), # Count peptides found in both replicates
    Percent_Peptides_In_One_Replicate = (Unique_Peptides_In_One_Replicate / Unique_Peptides) * 100, # Calculate % of peptides in one replicate
    Percent_Peptides_In_Both_Replicates = (Unique_Peptides_In_Both_Replicates / Unique_Peptides) * 100, # Calculate % of peptides in both replicates
    Abundance_Peptides_In_One_Replicate = sum(ifelse(Peptide_Presence == "one", Peptide_Abundance, 0), na.rm = TRUE), # Sum abundance of peptides in one replicate
    Abundance_Peptides_In_Both_Replicates = sum(ifelse(Peptide_Presence == "both", Peptide_Abundance, 0), na.rm = TRUE), # Sum abundance of peptides in both replicates
    Percent_Abundance_Peptides_In_One_Replicate = (Abundance_Peptides_In_One_Replicate / Total_Peptide_Abundance) * 100, # % abundance for one replicate peptides 
    Percent_Abundance_Peptides_In_Both_Replicates = (Abundance_Peptides_In_Both_Replicates / Total_Peptide_Abundance) * 100, # % abundance for both replicate peptides
    Overall_Unique_Peptides_Detected = first(Overall_Unique_Peptides_Detected), # fetch the overall unique peptides detected variable
    .groups = "drop") %>% # Drop grouping after summarizing
  mutate(Database = "NcbiRefSeq") # Add a column to indicate the database

# 3. Peptides per Sample per Dataset (NCBIrefseq)
Peptide_Summary_MockDB <- Peptide_Dat_MockDB %>% # Start with CostumDB peptides dataset
  filter(Peptide_Abundance > 0) %>% # Keep only rows with abundance > 0
  mutate(Overall_Unique_Peptides_Detected = n_distinct(Peptide_ID)) %>% # calculate the total unique peptides found in the whole dataset
  group_by(Date) %>% # Group by sample Date
  summarize(                                       
    Unique_Peptides = n_distinct(Peptide_ID), # Count unique peptides per sample
    Total_Peptide_Abundance = sum(Peptide_Abundance, na.rm = TRUE), # Sum total peptide abundance per sample
    Unique_Peptides_In_One_Replicate = sum(Peptide_Presence == "one"), # Count peptides found in one replicate
    Unique_Peptides_In_Both_Replicates = sum(Peptide_Presence == "both"), # Count peptides found in both replicates
    Percent_Peptides_In_One_Replicate = (Unique_Peptides_In_One_Replicate / Unique_Peptides) * 100, # Calculate % of peptides in one replicate
    Percent_Peptides_In_Both_Replicates = (Unique_Peptides_In_Both_Replicates / Unique_Peptides) * 100, # Calculate % of peptides in both replicates
    Abundance_Peptides_In_One_Replicate = sum(ifelse(Peptide_Presence == "one", Peptide_Abundance, 0), na.rm = TRUE), # Sum abundance of peptides in one replicate
    Abundance_Peptides_In_Both_Replicates = sum(ifelse(Peptide_Presence == "both", Peptide_Abundance, 0), na.rm = TRUE), # Sum abundance of peptides in both replicates
    Percent_Abundance_Peptides_In_One_Replicate = (Abundance_Peptides_In_One_Replicate / Total_Peptide_Abundance) * 100, # % abundance for one replicate peptides 
    Percent_Abundance_Peptides_In_Both_Replicates = (Abundance_Peptides_In_Both_Replicates / Total_Peptide_Abundance) * 100, # % abundance for both replicate peptides
    Overall_Unique_Peptides_Detected = first(Overall_Unique_Peptides_Detected), # fetch the overall unique peptides detected variable
    .groups = "drop") %>% # Drop grouping after summarizing
  mutate(Database = "MockDB") # Add a column to indicate the database

# 3. Merge datasets
Combined_Peptide_Summary <- bind_rows(Peptide_Summary_CostumDB, Peptide_Summary_UniProtDB, Peptide_Summary_NcbiRefSeq, Peptide_Summary_MockDB)  # Bind rows together

# ----------------------------Protein Dataset-----------------------------------
# 1. Proteins per Sample per Dataset (CostumDB)
Protein_Summary_CostumDB <- Protein_Dat_CostumDB %>% # Start with CostumDB proteins dataset
  filter(Protein_Abundance_IBAQ > 0) %>% # Keep only rows with abundance > 0
  mutate(Overall_Unique_Proteins_Detected = n_distinct(Protein_ID)) %>% # calculate the total unique proteins found in the whole dataset
  group_by(Date) %>% # Group by sample Date
  summarize(
    Unique_Proteins = n_distinct(Protein_ID), # Count unique proteins per sample
    Total_Protein_Abundance = sum(Protein_Abundance_IBAQ, na.rm = TRUE), # Sum total protein abundance per sample
    Unique_Proteins_In_One_Replicate = sum(Protein_Presence == "one"), # Count proteins found in one replicate
    Unique_Proteins_In_Both_Replicates = sum(Protein_Presence == "both"), # Count proteins found in both replicates
    Percent_Proteins_In_One_Replicate = (Unique_Proteins_In_One_Replicate / Unique_Proteins) * 100, # Calculate % of proteins in one replicate
    Percent_Proteins_In_Both_Replicates = (Unique_Proteins_In_Both_Replicates / Unique_Proteins) * 100, # Calculate % of proteins in both replicates
    Abundance_Proteins_In_One_Replicate = sum(ifelse(Protein_Presence == "one", Protein_Abundance_IBAQ, 0), na.rm = TRUE), # Sum abundance of proteins in one replicate
    Abundance_Proteins_In_Both_Replicates = sum(ifelse(Protein_Presence == "both", Protein_Abundance_IBAQ, 0), na.rm = TRUE), # Sum abundance of proteins in both replicates
    Percent_Abundance_Proteins_In_One_Replicate = (Abundance_Proteins_In_One_Replicate / Total_Protein_Abundance) * 100, # % abundance for one replicate proteins
    Percent_Abundance_Proteins_In_Both_Replicates = (Abundance_Proteins_In_Both_Replicates / Total_Protein_Abundance) * 100, # % abundance for both replicate proteins
    Overall_Unique_Proteins_Detected = first(Overall_Unique_Proteins_Detected), # fetch the overall unique peptides detected variable
    .groups = "drop") %>% # Drop grouping after summarizing
  mutate(Database = "CostumDB") # Add a column to indicate the database

# 2. Proteins per Sample per Dataset (UniProtDB)
Protein_Summary_UniProtDB <- Protein_Dat_UniProtDB %>% # Start with UniProtDB proteins dataset
  filter(Protein_Abundance_IBAQ > 0) %>% # Keep only rows with abundance > 0
  mutate(Overall_Unique_Proteins_Detected = n_distinct(Protein_ID)) %>% # calculate the total unique proteins found in the whole dataset
  group_by(Date) %>% # Group by sample Date
  summarize(
    Unique_Proteins = n_distinct(Protein_ID), # Count unique proteins per sample
    Total_Protein_Abundance = sum(Protein_Abundance_IBAQ, na.rm = TRUE), # Sum total protein abundance per sample
    Unique_Proteins_In_One_Replicate = sum(Protein_Presence == "one"), # Count proteins found in one replicate
    Unique_Proteins_In_Both_Replicates = sum(Protein_Presence == "both"), # Count proteins found in both replicates
    Percent_Proteins_In_One_Replicate = (Unique_Proteins_In_One_Replicate / Unique_Proteins) * 100, # Calculate % of proteins in one replicate
    Percent_Proteins_In_Both_Replicates = (Unique_Proteins_In_Both_Replicates / Unique_Proteins) * 100, # Calculate % of proteins in both replicates
    Abundance_Proteins_In_One_Replicate = sum(ifelse(Protein_Presence == "one", Protein_Abundance_IBAQ, 0), na.rm = TRUE), # Sum abundance of proteins in one replicate
    Abundance_Proteins_In_Both_Replicates = sum(ifelse(Protein_Presence == "both", Protein_Abundance_IBAQ, 0), na.rm = TRUE), # Sum abundance of proteins in both replicates
    Percent_Abundance_Proteins_In_One_Replicate = (Abundance_Proteins_In_One_Replicate / Total_Protein_Abundance) * 100, # % abundance for one replicate proteins
    Percent_Abundance_Proteins_In_Both_Replicates = (Abundance_Proteins_In_Both_Replicates / Total_Protein_Abundance) * 100, # % abundance for both replicate proteins
    Overall_Unique_Proteins_Detected = first(Overall_Unique_Proteins_Detected), # fetch the overall unique peptides detected variable
    .groups = "drop") %>% # Drop grouping after summarizing
  mutate(Database = "UniProtDB") # Add a column to indicate the database

# 2. Proteins per Sample per Dataset (NCBIREfSeq)
Protein_Summary_NcbiRefSeq <- Protein_Dat_NcbiRefSeq %>% # Start with UniProtDB proteins dataset
  filter(Protein_Abundance_IBAQ > 0) %>% # Keep only rows with abundance > 0
  mutate(Overall_Unique_Proteins_Detected = n_distinct(Protein_ID)) %>% # calculate the total unique proteins found in the whole dataset
  group_by(Date) %>% # Group by sample Date
  summarize(
    Unique_Proteins = n_distinct(Protein_ID), # Count unique proteins per sample
    Total_Protein_Abundance = sum(Protein_Abundance_IBAQ, na.rm = TRUE), # Sum total protein abundance per sample
    Unique_Proteins_In_One_Replicate = sum(Protein_Presence == "one"), # Count proteins found in one replicate
    Unique_Proteins_In_Both_Replicates = sum(Protein_Presence == "both"), # Count proteins found in both replicates
    Percent_Proteins_In_One_Replicate = (Unique_Proteins_In_One_Replicate / Unique_Proteins) * 100, # Calculate % of proteins in one replicate
    Percent_Proteins_In_Both_Replicates = (Unique_Proteins_In_Both_Replicates / Unique_Proteins) * 100, # Calculate % of proteins in both replicates
    Abundance_Proteins_In_One_Replicate = sum(ifelse(Protein_Presence == "one", Protein_Abundance_IBAQ, 0), na.rm = TRUE), # Sum abundance of proteins in one replicate
    Abundance_Proteins_In_Both_Replicates = sum(ifelse(Protein_Presence == "both", Protein_Abundance_IBAQ, 0), na.rm = TRUE), # Sum abundance of proteins in both replicates
    Percent_Abundance_Proteins_In_One_Replicate = (Abundance_Proteins_In_One_Replicate / Total_Protein_Abundance) * 100, # % abundance for one replicate proteins
    Percent_Abundance_Proteins_In_Both_Replicates = (Abundance_Proteins_In_Both_Replicates / Total_Protein_Abundance) * 100, # % abundance for both replicate proteins
    Overall_Unique_Proteins_Detected = first(Overall_Unique_Proteins_Detected), # fetch the overall unique peptides detected variable
    .groups = "drop") %>% # Drop grouping after summarizing
  mutate(Database = "NcbiRefSeq") # Add a column to indicate the database

# 2. Proteins per Sample per Dataset (NCBIREfSeq)
Protein_Summary_MockDB <- Protein_Dat_MockDB %>% # Start with UniProtDB proteins dataset
  filter(Protein_Abundance_IBAQ > 0) %>% # Keep only rows with abundance > 0
  mutate(Overall_Unique_Proteins_Detected = n_distinct(Protein_ID)) %>% # calculate the total unique proteins found in the whole dataset
  group_by(Date) %>% # Group by sample Date
  summarize(
    Unique_Proteins = n_distinct(Protein_ID), # Count unique proteins per sample
    Total_Protein_Abundance = sum(Protein_Abundance_IBAQ, na.rm = TRUE), # Sum total protein abundance per sample
    Unique_Proteins_In_One_Replicate = sum(Protein_Presence == "one"), # Count proteins found in one replicate
    Unique_Proteins_In_Both_Replicates = sum(Protein_Presence == "both"), # Count proteins found in both replicates
    Percent_Proteins_In_One_Replicate = (Unique_Proteins_In_One_Replicate / Unique_Proteins) * 100, # Calculate % of proteins in one replicate
    Percent_Proteins_In_Both_Replicates = (Unique_Proteins_In_Both_Replicates / Unique_Proteins) * 100, # Calculate % of proteins in both replicates
    Abundance_Proteins_In_One_Replicate = sum(ifelse(Protein_Presence == "one", Protein_Abundance_IBAQ, 0), na.rm = TRUE), # Sum abundance of proteins in one replicate
    Abundance_Proteins_In_Both_Replicates = sum(ifelse(Protein_Presence == "both", Protein_Abundance_IBAQ, 0), na.rm = TRUE), # Sum abundance of proteins in both replicates
    Percent_Abundance_Proteins_In_One_Replicate = (Abundance_Proteins_In_One_Replicate / Total_Protein_Abundance) * 100, # % abundance for one replicate proteins
    Percent_Abundance_Proteins_In_Both_Replicates = (Abundance_Proteins_In_Both_Replicates / Total_Protein_Abundance) * 100, # % abundance for both replicate proteins
    Overall_Unique_Proteins_Detected = first(Overall_Unique_Proteins_Detected), # fetch the overall unique peptides detected variable
    .groups = "drop") %>% # Drop grouping after summarizing
  mutate(Database = "MockDB") # Add a column to indicate the database

# 3. Merge protein datasets
Combined_Protein_Summary <- bind_rows(Protein_Summary_CostumDB, Protein_Summary_UniProtDB, Protein_Summary_NcbiRefSeq, Protein_Summary_MockDB) # Bind rows together to create one dataset


################################################################################
#                              EXPORT DATASET                                  #
################################################################################
# Merge the Datasets
DatabaseComparison_Summary <- Combined_Peptide_Summary %>%
  left_join(Combined_Protein_Summary, by = c("Date", "Database"))

# Write the data to csv
write_csv(DatabaseComparison_Summary, "DatabaseComparison_Summary_FourDatabases.csv") # Export combined dataset to CSV

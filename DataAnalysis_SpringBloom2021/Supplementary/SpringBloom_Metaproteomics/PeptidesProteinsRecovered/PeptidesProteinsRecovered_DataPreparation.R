################################################################################
#                                                                              #
#                             Peptides/Proteins Recovered                      #
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
library(readxl)

################################################################################
#                              SET WORKING DIRECTORY                           #
################################################################################
# Adjust the path as per your project location
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Method_Paper_Analysis/PeptidesProteinsRecovered/")


################################################################################
#                              LOAD DATASETS                                   #
################################################################################
diann <- read_delim("../../Data/Proteomics_Data/Raw_Data/72Samples_p34548_SpringBloom2021_95MagSpeciesDB_15JUL24_20240726_Database_Run07JAN25/diann-report.tsv")

Protein_Dat <- read_csv("../../Datasets/72Samples_CostumDB/Integrated_Protein_Table.csv", col_select = c(Protein_ID, Species_ID, KEGG_ID, Phylum, Protein_Abundance_IBAQ, Date))
Peptide_Dat <- read_csv("../../Datasets/72Samples_CostumDB/Integrated_Peptide_Table.csv", col_select = c(Protein_ID, Species_ID, Peptide_Sequence, Peptide_ID, Protein_Sequence))


################################################################################
#                              PEPARE DATASETS                                 #
################################################################################
############################################################

# Define the four stages of filtering
f1 <- diann %>% filter(Q.Value       < 0.01)
f2 <- diann %>% filter(Q.Value       < 0.01,
                       Lib.Q.Value   < 0.01,
                       PG.Q.Value    < 0.01)
f3 <- diann %>% filter(Q.Value       < 0.01,
                       Lib.Q.Value   < 0.01,
                       PG.Q.Value    < 0.01) %>%
  group_by(Protein.Group) %>%
  filter(n_distinct(Stripped.Sequence) >= 2) %>%
  ungroup()

stages <- list(
  Raw      = diann,
  `Q<1%`   = f1,
  `All<1%` = f2,
  `2+ pep` = f3
)

# Summarize counts at each stage
summary_tbl <- imap_dfr(stages, ~ {
  df <- .x
  # how many unique peptides & which of them are proteotypic?
  pep_ids  <- df %>% distinct(Stripped.Sequence, Proteotypic)
  # how many unique protein-groups & which have a proteotypic anchor?
  prot_ids <- df %>% distinct(Protein.Group, Stripped.Sequence, Proteotypic)
  
  tibble(
    Stage     = .y,
    Peptides  = n_distinct(df$Stripped.Sequence),
    Proteins  = n_distinct(df$Protein.Group),
    ProteotypicPeptides   = sum(pep_ids$Proteotypic == 1),
    ProteotypicProteins = prot_ids %>%
      filter(Proteotypic == 1) %>%
      pull(Protein.Group) %>%
      unique() %>%
      length()
  )
}) %>%
  mutate(Stage = factor(Stage, levels = c("Raw","Q<1%","All<1%","2+ pep")))


# prepare the data for the length overview
dat <- Protein_Dat %>%
  filter(Phylum != "Y-FGCZCont", Protein_Abundance_IBAQ > 0) %>%
  select(-Protein_Abundance_IBAQ, -Date) %>%
  distinct() %>%
  left_join(Peptide_Dat, by = c("Protein_ID", "Species_ID")) %>%
  mutate(Peptide_Length = nchar(Peptide_Sequence),
         Protein_Length = nchar(Protein_Sequence))

################################################################################
#                              EXPORT DATASET                                  #
################################################################################
# Write the data to csv
write_csv(summary_tbl, "PeptideProteinRecovery_Summary.csv")
write_csv(dat, "PeptideProteinDistribution.csv")


test <- Protein_Dat %>%
  filter(Protein_Abundance_IBAQ > 0) %>%
  group_by(Date) %>%
  summarize(n_distinct(Protein_ID))

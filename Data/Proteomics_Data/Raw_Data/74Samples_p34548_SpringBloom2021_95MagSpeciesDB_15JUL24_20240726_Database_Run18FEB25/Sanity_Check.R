rm(list = ls())

library(readxl)
library(dplyr)
library(ggplot2)

#set working directory
setwd("C:/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Proteomics_Data/Raw_Data/74Samples_p34548_SpringBloom2021_95MagSpeciesDB_15JUL24_20240726_Database_Run18FEB25/")


##################################################

# read the data
prolfqua <- read_excel("proteinAbundances_.xlsx", sheet = 4)

# summary all peptides and proteins
prolfqua_allproteins <- prolfqua %>% 
  select(protein_Id, nrPeptides) %>%
  filter(nrPeptides >= 1) %>%
  distinct() %>%
  summarise(Unique_Peptides = sum(nrPeptides),
            Unique_Proteins = n_distinct(protein_Id))

# summary more than two peptides per protein
prolfqua_more2peptides <- prolfqua %>% 
  select(protein_Id, nrPeptides) %>%
  filter(nrPeptides >= 2) %>%
  distinct() %>%
  summarise(Unique_Peptides = sum(nrPeptides),
            Unique_Proteins = n_distinct(protein_Id))

############################################################
diann <- read.delim("diann-report.tsv")

# original
diann_summary <- diann %>%
  select(Protein.Group, Stripped.Sequence) %>%
  summarise(Unique_Peptides = n_distinct(Stripped.Sequence),
            Unique_Proteins = n_distinct(Protein.Group))

# filter at 0.01 Q.Value
diann_summary_1FDRPeptides <- diann %>%
  select(Protein.Group, Stripped.Sequence, Q.Value) %>%
  filter(Q.Value < 0.01) %>%
  summarise(Unique_Peptides = n_distinct(Stripped.Sequence),
            Unique_Proteins = n_distinct(Protein.Group))

# filter at Lib.PG.Q.Value < 0.01, PG.Q.Value < 0.01 and Q.Value < 0.01
diann_summary_1FDRPeptides_1FDRProteins <- diann %>%
  select(Protein.Group, Stripped.Sequence, Lib.PG.Q.Value, PG.Q.Value) %>%
  filter(Lib.PG.Q.Value < 0.01 & PG.Q.Value < 0.01) %>%
  summarise(Unique_Peptides = n_distinct(Stripped.Sequence),
            Unique_Proteins = n_distinct(Protein.Group))

# additionally filter proteins with less than two peptides
diann_summary_1FDRPeptides_1FDRProteins_2Peptides <- diann %>%
  select(Protein.Group, Stripped.Sequence, Lib.PG.Q.Value, PG.Q.Value, Q.Value) %>%
  filter(Lib.PG.Q.Value < 0.01 & PG.Q.Value < 0.01 & Q.Value < 0.01) %>%
  group_by(Protein.Group) %>%
  filter(n_distinct(Stripped.Sequence) >= 2) %>%
  ungroup() %>%
  distinct() %>%
  summarise(Unique_Peptides = n_distinct(Stripped.Sequence),
            Unique_Proteins = n_distinct(Protein.Group))

# export sanity check summary:
write.csv(diann_summary, "Diann_Raw.csv", row.names = FALSE)
write.csv(diann_summary_1FDRPeptides_1FDRProteins_2Peptides, "Diann_Final.csv", row.names = FALSE)
write.csv(prolfqua_allproteins, "Prolefqua_Raw.csv",row.names = FALSE)
write.csv(prolfqua_more2peptides, "Prolefqua_Final.csv",row.names = FALSE)

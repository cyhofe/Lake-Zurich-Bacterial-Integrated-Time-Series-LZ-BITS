rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(stringr)

setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Figure_02/ProteomeRecoveryVsSpeciesAbundance/")

# --- Load data ---
Protein_Dat <- read_csv(
  "../../../Datasets/72Samples_CostumDB/Integrated_Protein_Table.csv",
  col_select = c(
    Date, Protein_ID, Species_ID, Protein_Abundance_IBAQ,
    Phylum, Class, Order, Family, Genus, Species,
    Estimated_Genome_Size, Protein_Group
  )
)

Species_Dat <- read_csv("../../../Data/Species_Data/Final_Data/Species_Table.csv") %>%
  select(Species_ID, Date, Species_Abundance, Phylum, Class, Order, Family, Genus, Species, Estimated_Genome_Size)

Protein_Dat <- Protein_Dat %>% mutate(Date = as.Date(Date))
Species_Dat <- Species_Dat %>% mutate(Date = as.Date(Date))

# =============================================================================
# 1) Remove contaminants
# =============================================================================
Protein_Dat <- Protein_Dat %>%
  filter(Phylum != "Y-FGCZCont")

# =============================================================================
# 2) Median-normalize protein abundances (per Date)
# =============================================================================
sample_medians <- Protein_Dat %>%
  filter(!is.na(Protein_Abundance_IBAQ), Protein_Abundance_IBAQ > 0) %>%
  group_by(Date) %>%
  summarise(med_iBAQ = median(Protein_Abundance_IBAQ, na.rm = TRUE), .groups = "drop")

global_med <- median(sample_medians$med_iBAQ, na.rm = TRUE)

Protein_Dat <- Protein_Dat %>%
  left_join(sample_medians, by = "Date") %>%
  mutate(Protein_Abundance_IBAQ_Norm = Protein_Abundance_IBAQ / med_iBAQ * global_med) %>%
  select(-med_iBAQ)

# =============================================================================
# 3) Compute species-specific proteins (based on Protein_Group taxonomy consistency)
# =============================================================================
protein_cleaned <- Protein_Dat %>%
  filter(!is.na(Protein_Abundance_IBAQ_Norm), Protein_Abundance_IBAQ_Norm > 0) %>%
  select(Protein_ID, Protein_Group) %>%
  distinct()

taxonomy_matches <- protein_cleaned %>%
  separate_rows(Protein_Group, sep = ";") %>%
  rename(Rep_Protein_ID = Protein_ID,
         Group_Member_ID = Protein_Group) %>%
  mutate(
    Rep_Species_ID  = str_extract(Rep_Protein_ID, "^[^-]+"),
    Memb_Species_ID = str_extract(Group_Member_ID, "^[^-]+")
  ) %>%
  left_join(
    Species_Dat %>% distinct(Species_ID, Phylum, Class, Order, Family, Genus, Species),
    by = c("Rep_Species_ID" = "Species_ID")
  ) %>%
  rename_with(~ paste0(.x, "_Rep"), c(Phylum, Class, Order, Family, Genus, Species)) %>%
  left_join(
    Species_Dat %>% distinct(Species_ID, Phylum, Class, Order, Family, Genus, Species),
    by = c("Memb_Species_ID" = "Species_ID")
  ) %>%
  rename_with(~ paste0(.x, "_Memb"), c(Phylum, Class, Order, Family, Genus, Species)) %>%
  mutate(
    Species_Match = Species_Rep == Species_Memb,
    Genus_Match   = Genus_Rep   == Genus_Memb,
    Family_Match  = Family_Rep  == Family_Memb,
    Order_Match   = Order_Rep   == Order_Memb,
    Class_Match   = Class_Rep   == Class_Memb,
    Phylum_Match  = Phylum_Rep  == Phylum_Memb
  ) %>%
  mutate(
    Resolution = case_when(
      Species_Match ~ "Species",
      Genus_Match   ~ "Genus",
      Family_Match  ~ "Family",
      Order_Match   ~ "Order",
      Class_Match   ~ "Class",
      Phylum_Match  ~ "Phylum",
      TRUE          ~ "NoMatch"
    )
  )

resolution_summary <- taxonomy_matches %>%
  group_by(Rep_Protein_ID) %>%
  summarise(
    Species_Count = sum(Resolution == "Species"),
    Genus_Count   = sum(Resolution == "Genus"),
    Family_Count  = sum(Resolution == "Family"),
    Order_Count   = sum(Resolution == "Order"),
    Class_Count   = sum(Resolution == "Class"),
    Phylum_Count  = sum(Resolution == "Phylum"),
    NoMatch_Count = sum(Resolution == "NoMatch"),
    .groups = "drop"
  )

species_specific_proteins <- resolution_summary %>%
  filter(
    Species_Count > 0,
    Genus_Count   == 0,
    Family_Count  == 0,
    Order_Count   == 0,
    Class_Count   == 0,
    Phylum_Count  == 0,
    NoMatch_Count == 0
  ) %>%
  pull(Rep_Protein_ID) %>%
  unique()

# =============================================================================
# 4) Keep only species-specific proteins
# =============================================================================
Protein_Dat_SS <- Protein_Dat %>%
  filter(Protein_ID %in% species_specific_proteins)

# =============================================================================
# 5) Remove entries not detected (Norm == 0) and where species abundance == 0
#    (explicitly join metagenome abundance per Species_ID+Date)
# =============================================================================
Protein_Dat_SS <- Protein_Dat_SS %>%
  left_join(
    Species_Dat %>% select(Species_ID, Date, Species_Abundance),
    by = c("Species_ID", "Date")
  ) %>%
  filter(
    !is.na(Protein_Abundance_IBAQ_Norm),
    Protein_Abundance_IBAQ_Norm > 0,
    !is.na(Species_Abundance),
    Species_Abundance > 0
  )

# =============================================================================
# 6) For each species: select date where species abundance is maximal (ties -> earliest)
# =============================================================================
species_peak <- Protein_Dat_SS %>%
  distinct(Species_ID, Date, Species_Abundance) %>%
  group_by(Species_ID) %>%
  filter(Species_Abundance == max(Species_Abundance, na.rm = TRUE)) %>%
  arrange(Date) %>%
  slice(1) %>%
  ungroup() %>%
  rename(Peak_Date = Date,
         Species_Abundance_at_Peak = Species_Abundance)

# =============================================================================
# 7) Summarize ON THAT DATE: total protein signal + number of proteins detected
# =============================================================================
species_summary_peak <- Protein_Dat_SS %>%
  inner_join(species_peak, by = "Species_ID") %>%
  filter(Date == Peak_Date) %>%
  group_by(Species_ID) %>%
  summarise(
    Peak_Date = first(Peak_Date),
    Species_Abundance_at_Peak = first(Species_Abundance_at_Peak),
    Proteome_Abundance_on_PeakDate = sum(Protein_Abundance_IBAQ_Norm, na.rm = TRUE),
    Detected_Proteins_at_Peak = n_distinct(Protein_ID),
    .groups = "drop"
  ) %>%
  left_join(
    Species_Dat %>%
      distinct(Species_ID, Phylum, Class, Order, Family, Genus, Species, Estimated_Genome_Size),
    by = "Species_ID"
  )

# --- Export ---
write.csv(species_summary_peak, "species_summary_peak_only.csv", row.names = FALSE)

# --- Quick sanity output ---
cat("Species-specific proteins used:", length(species_specific_proteins), "\n")
cat("Species in final summary:      ", n_distinct(species_summary_peak$Species_ID), "\n")
cat("Rows in final summary:         ", nrow(species_summary_peak), "\n")

################################################################################
#                                                                              #
#                             Proteome Specificity                             #
#                                                                              #
#   Author  : Cyrill Hofer                                                     #
#   Date    : 09.04.2025                                                       #
#   Purpose : Data Preparation (with Date-wise median-centering normalization)  #
#                                                                              #
################################################################################
################################################################################
# Proteome Specificity — FIXED & ROCK-SOLID protein-group resolution
# - Repairs missing columns / wrong grouping
# - Computes LCA-style resolution from unfolded Protein_Group members
# - Produces reviewer-safe Species_Specific_ProteinGroup
################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Figure_02/ProteomeSpecificity/")

# ------------------------------ LOAD DATASETS ---------------------------------
Protein_Dat <- read_csv(
  "../../../Datasets/72Samples_CostumDB/Integrated_Protein_Table.csv",
  col_select = c(Date, Protein_ID, Species_ID, Protein_Group, Protein_Abundance_IBAQ, Phylum),
  show_col_types = FALSE
)

Species_Dat <- read_csv(
  "../../../Data/Species_Data/Final_Data/Species_Table.csv",
  col_select = c(Species_ID, Phylum, Class, Order, Family, Genus, Species),
  show_col_types = FALSE
) %>% distinct()

# remove contaminants
Protein_Dat <- Protein_Dat %>%
  filter(Phylum != "Y-FGCZCont")

# ---------------- Date-wise median-centering normalization ---------------------
sample_medians <- Protein_Dat %>%
  filter(Protein_Abundance_IBAQ > 0) %>%
  group_by(Date) %>%
  summarise(med_iBAQ = median(Protein_Abundance_IBAQ, na.rm = TRUE), .groups = "drop")

global_med <- median(sample_medians$med_iBAQ, na.rm = TRUE)

Protein_Dat <- Protein_Dat %>%
  left_join(sample_medians, by = "Date") %>%
  mutate(Protein_Abundance_IBAQ_Norm = Protein_Abundance_IBAQ / med_iBAQ * global_med) %>%
  select(-med_iBAQ)

# ------------------------------ PREPARE PROTEIN TABLE -------------------------
protein_dat_prepared <- Protein_Dat %>%
  filter(
    Phylum != "Y-FGCZCont",
    Protein_Abundance_IBAQ > 0
  ) %>%
  group_by(Protein_ID, Species_ID, Protein_Group) %>%
  summarise(
    Total_Abundance_Per_Protein_Raw = sum(Protein_Abundance_IBAQ, na.rm = TRUE),
    Total_Abundance_Per_Protein     = sum(Protein_Abundance_IBAQ_Norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Protein_Group_Size = str_count(Protein_Group, ";") + 1)

# ==============================================================================
# FIXED: Build a Protein_ID-level LUT with LCA-style resolution
# ------------------------------------------------------------------------------
# Definitions:
#   - Representative = Protein_ID (your curated representative protein)
#   - Group members  = identifiers in Protein_Group (split by ';')
#   - Resolution     = finest rank where ALL group members share the same label
#                      (Species -> ... -> Phylum), else No_Match.
#
# Conservative rule (as in your DIA-NN script):
#   - if ANY group member cannot be mapped to taxonomy -> No_Match
# ==============================================================================

# 1) Expand members and join taxonomy for representative + members
taxonomy_pairs <- protein_dat_prepared %>%
  mutate(
    Species_ID_Representative = str_extract(Protein_ID, "^[^-]+")
  ) %>%
  separate_rows(Protein_Group, sep = ";") %>%
  mutate(
    Protein_Group_Member   = str_trim(Protein_Group),
    Species_ID_GroupMember = str_extract(Protein_Group_Member, "^[^-]+")
  ) %>%
  # join representative taxonomy
  left_join(Species_Dat, by = c("Species_ID_Representative" = "Species_ID")) %>%
  rename_with(~ paste0(.x, "_Rep"), c(Phylum, Class, Order, Family, Genus, Species)) %>%
  # join group-member taxonomy
  left_join(Species_Dat, by = c("Species_ID_GroupMember" = "Species_ID")) %>%
  rename_with(~ paste0(.x, "_Grp"), c(Phylum, Class, Order, Family, Genus, Species))

# 2) Summarise per Protein_ID (your curated “protein group” unit)
protein_group_lut <- taxonomy_pairs %>%
  group_by(Protein_ID, Species_ID_Representative) %>%
  summarise(
    # carry representative taxonomy (stable per Protein_ID)
    Phylum_Rep  = first(Phylum_Rep),
    Class_Rep   = first(Class_Rep),
    Order_Rep   = first(Order_Rep),
    Family_Rep  = first(Family_Rep),
    Genus_Rep   = first(Genus_Rep),
    Species_Rep = first(Species_Rep),
    
    # how many distinct taxa among GROUP MEMBERS (identity-aware)
    n_species = n_distinct(Species_Grp, na.rm = TRUE),
    n_genus   = n_distinct(Genus_Grp,   na.rm = TRUE),
    n_family  = n_distinct(Family_Grp,  na.rm = TRUE),
    n_order   = n_distinct(Order_Grp,   na.rm = TRUE),
    n_class   = n_distinct(Class_Grp,   na.rm = TRUE),
    n_phylum  = n_distinct(Phylum_Grp,  na.rm = TRUE),
    
    # conservative: if any member is unmapped at phylum => treat group as No_Match
    any_missing_taxonomy = any(is.na(Phylum_Grp)),
    
    .groups = "drop"
  ) %>%
  mutate(
    ProteinGroup_Resolution = case_when(
      any_missing_taxonomy ~ "No_Match",
      n_species == 1 ~ "Species_Resolution",
      n_genus   == 1 ~ "Genus_Resolution",
      n_family  == 1 ~ "Family_Resolution",
      n_order   == 1 ~ "Order_Resolution",
      n_class   == 1 ~ "Class_Resolution",
      n_phylum  == 1 ~ "Phylum_Resolution",
      TRUE           ~ "No_Match"
    ),
    Species_Specific_ProteinGroup = (ProteinGroup_Resolution == "Species_Resolution")
  )

# 3) Merge LUT back to your curated protein table
dat <- protein_dat_prepared %>%
  left_join(
    protein_group_lut %>%
      rename(Species_ID = Species_ID_Representative),
    by = c("Protein_ID", "Species_ID")
  )

# ------------------------------ EXPORT ----------------------------------------
write_csv(dat, "ProteinGroupResolution.csv")

# ==============================================================================
# Summary table: how many proteins at each specificity level
# Assumes your final merged table is called `dat` and contains:
#   - Protein_ID
#   - ProteinGroup_Resolution
#   - Species_Specific_ProteinGroup (optional)
#   - Total_Abundance_Per_Protein (optional; for weighted summaries)
# ==============================================================================

library(dplyr)
library(tidyr)

# ---------- 1) Simple counts (unique proteins per resolution class) ----------
summary_specificity_counts <- dat %>%
  distinct(Protein_ID, ProteinGroup_Resolution) %>%
  count(ProteinGroup_Resolution, name = "n_proteins") %>%
  mutate(
    frac_proteins = n_proteins / sum(n_proteins)
  ) %>%
  arrange(desc(n_proteins))

summary_specificity_counts


# ---------- 2) Same, but with explicit ordered factor (nice for plotting/reporting) ----------
resolution_levels <- c(
  "Species_Resolution", "Genus_Resolution", "Family_Resolution",
  "Order_Resolution", "Class_Resolution", "Phylum_Resolution", "No_Match"
)

summary_specificity_counts_ordered <- dat %>%
  distinct(Protein_ID, ProteinGroup_Resolution) %>%
  mutate(ProteinGroup_Resolution = factor(ProteinGroup_Resolution, levels = resolution_levels)) %>%
  count(ProteinGroup_Resolution, name = "n_proteins") %>%
  mutate(frac_proteins = n_proteins / sum(n_proteins)) %>%
  arrange(ProteinGroup_Resolution)

summary_specificity_counts_ordered


# ---------- 3) Optional: abundance-weighted summary (if you want “proteome share” per class) ----------
# This answers: "what fraction of total protein abundance is assigned to each specificity level?"
# (Uses Total_Abundance_Per_Protein; change to *_Raw if you prefer.)
if ("Total_Abundance_Per_Protein" %in% colnames(dat)) {
  
  summary_specificity_abundance <- dat %>%
    distinct(Protein_ID, ProteinGroup_Resolution, Total_Abundance_Per_Protein) %>%
    group_by(ProteinGroup_Resolution) %>%
    summarise(
      n_proteins = n(),
      sum_abundance = sum(Total_Abundance_Per_Protein, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      frac_proteins = n_proteins / sum(n_proteins),
      frac_abundance = sum_abundance / sum(sum_abundance)
    ) %>%
    arrange(desc(sum_abundance))
  
  summary_specificity_abundance
}


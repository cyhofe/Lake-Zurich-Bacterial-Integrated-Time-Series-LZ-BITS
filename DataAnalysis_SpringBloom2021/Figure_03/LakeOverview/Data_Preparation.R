############################################################
# Lake Overview Data Preparation
#
# Purpose
#   Prepare analysis-ready tables for joint metagenomic and metaproteomic
#   analyses while preserving the FULL species universe.
#
# Key design decisions 
#   1) Keep ALL species (MAGs) in species_meta (no silent drops).
#   2) Restrict metagenomic abundance to proteomics sampling dates to align
#      genomes and proteomes in time.
#   3) Encode missing metagenomic observations on proteomics dates as explicit 0
#      (so "absent" ≠ "missing row").
#   4) Keep only species-specific protein groups (conservative taxonomy-based rule).
#   5) Use ONE grouping variable for plots and downstream summaries:
#        Group_For_Plots
#      A species is assigned its trophic group ONLY if:
#        - TrophicIndex is available (proxy for >=80% completeness)
#        - Proteome_Recovered == TRUE
#        - Abundance_In_Proteomics_Dates == TRUE
#      Otherwise it is labelled: "Filtered"
#
# Output files
#   - Species_Meta_Table.csv
#   - Species_Abundance_Table.csv
#   - Species_Protein_Table.csv
#   - Species_Gene_Table.csv
############################################################

# ============================
# 1) Clear environment
# ============================
rm(list = ls())
gc()

# ============================
# 2) Libraries
# ============================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

# ============================
# 3) Working directory
# ============================
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Figure_03/LakeOverview/")

# ============================
# 4) Load input data
# ============================
species_dat <- read_csv("../../../Data/Species_Data/Final_Data/Species_Table.csv")
trophic_dat <- read_csv("../TrophicIndex/LakeZurich_RetainedSpecies_TrophicIndex_TrophicGroup_EST.csv")
kegg_dat    <- read_csv("../../../Data/KEGG_Category_Data/Final_Data/KEGG_Category_Table.csv")
gene_dat    <- read_csv("../../../Data/Gene_Data/Final_Data/Gene_Table.csv")
protein_dat <- read_csv("../../../Datasets/72Samples_CostumDB/Integrated_Protein_Table.csv")

############################################################
# A) Canonical species metadata (ALL species; no silent drops)
############################################################

# Unique taxonomy per Species_ID
species_tax <- species_dat %>%
  select(Species_ID, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  distinct()

# Trophic information
# Biological interpretation:
#   trophic indices were computed only for MAGs meeting the >=80% completeness
#   criterion in the trophic index workflow. Therefore, missing TrophicIndex is
#   interpreted as <80% completeness (or excluded from that derivation).
species_trophic <- trophic_dat %>%
  select(Species_ID, TrophicIndex, TrophicGroup) %>%
  distinct() %>%
  mutate(
    TrophicGroup = ifelse(is.na(TrophicGroup) | TrophicGroup == "", "not_defined", TrophicGroup)
  )

# Canonical species-level metadata table
species_meta <- species_tax %>%
  left_join(species_trophic, by = "Species_ID") %>%
  mutate(
    TrophicInfo_Available = !is.na(TrophicIndex),
    TrophicGroup = ifelse(is.na(TrophicGroup) | TrophicGroup == "", "not_defined", TrophicGroup)
  )

############################################################
# B) Date universe for JOINT analyses = proteomics sampling dates
############################################################
proteomics_dates <- sort(unique(protein_dat$Date))

############################################################
# C) Species abundance (Date-resolved, explicit zeros, proteomics dates only)
############################################################

# Build complete Species_ID × Date grid for proteomics dates.
# Missing abundance observations are encoded as 0 (explicit absence).
species_abundance <- tidyr::expand_grid(
  Species_ID = unique(species_meta$Species_ID),
  Date       = proteomics_dates
) %>%
  left_join(
    species_dat %>%
      filter(Date %in% proteomics_dates) %>%
      select(Species_ID, Date, Species_Abundance) %>%
      group_by(Species_ID, Date) %>%
      summarise(Species_Abundance = sum(Species_Abundance, na.rm = TRUE), .groups = "drop"),
    by = c("Species_ID", "Date")
  ) %>%
  mutate(Species_Abundance = ifelse(is.na(Species_Abundance), 0, Species_Abundance))

############################################################
# D) Protein table (Date-resolved; contaminants removed; median-centred; species-specific only)
############################################################

# Taxonomy lookup for protein-group specificity testing
species_taxonomy <- species_dat %>%
  select(Species_ID, Phylum, Class, Order, Family, Genus, Species) %>%
  distinct()

# ---- D1) Remove contaminants + align dates ----
Protein_Dat <- protein_dat %>%
  select(Date, Species_ID, Protein_ID, Protein_Group, Protein_Abundance_IBAQ, Phylum) %>%
  filter(Phylum != "Y-FGCZCont") %>%
  filter(Date %in% proteomics_dates)

# ---- D2) Date-wise median-centering normalisation (EXACT logic) ----
sample_medians <- Protein_Dat %>%
  filter(Protein_Abundance_IBAQ > 0) %>%
  group_by(Date) %>%
  summarise(med_iBAQ = median(Protein_Abundance_IBAQ, na.rm = TRUE), .groups = "drop")

global_med <- median(sample_medians$med_iBAQ, na.rm = TRUE)

Protein_Dat <- Protein_Dat %>%
  left_join(sample_medians, by = "Date") %>%
  mutate(Protein_Abundance_IBAQ_Norm = Protein_Abundance_IBAQ / med_iBAQ * global_med) %>%
  select(-med_iBAQ)

# ---- D3) Date-resolved protein abundance table ----
protein_date_table <- Protein_Dat %>%
  filter(Protein_Abundance_IBAQ > 0) %>%
  group_by(Date, Species_ID, Protein_ID, Protein_Group) %>%
  summarise(
    Protein_Abundance_IBAQ_Raw  = sum(Protein_Abundance_IBAQ, na.rm = TRUE),
    Protein_Abundance_IBAQ_Norm = sum(Protein_Abundance_IBAQ_Norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Protein_Group_Size = str_count(Protein_Group, ";") + 1)

# ---- D4) Species-specific protein groups (conservative taxonomy-based rule) ----
# Species-specific = all group members resolve to exactly one Species label.
# If any member is unmapped => not species-specific.
protein_group_specificity <- protein_date_table %>%
  distinct(Protein_ID, Protein_Group) %>%
  mutate(Species_ID_Representative = str_extract(Protein_ID, "^[^-]+")) %>%
  separate_rows(Protein_Group, sep = ";") %>%
  mutate(
    Protein_Group_Member   = str_trim(Protein_Group),
    Species_ID_GroupMember = str_extract(Protein_Group_Member, "^[^-]+")
  ) %>%
  left_join(
    species_taxonomy %>% select(Species_ID, Species),
    by = c("Species_ID_GroupMember" = "Species_ID")
  ) %>%
  group_by(Protein_ID, Species_ID_Representative) %>%
  summarise(
    any_missing_taxonomy = any(is.na(Species)),
    n_species_labels     = n_distinct(Species, na.rm = TRUE),
    Species_Specific_ProteinGroup = (!any_missing_taxonomy & n_species_labels == 1),
    .groups = "drop"
  )

# ---- D5) Final protein table: species-specific only ----
protein_table_final <- protein_date_table %>%
  left_join(
    protein_group_specificity %>% rename(Species_ID = Species_ID_Representative),
    by = c("Protein_ID", "Species_ID")
  ) %>%
  filter(Species_Specific_ProteinGroup %in% TRUE)

############################################################
# E) QC flags (computed once; stored in species_meta; no silent drops)
############################################################

# Proteome recovered per species = any detected species-specific protein signal
proteome_by_species <- protein_table_final %>%
  group_by(Species_ID) %>%
  summarise(
    Proteome_Total_iBAQ_Norm = sum(Protein_Abundance_IBAQ_Norm, na.rm = TRUE),
    Proteome_Recovered       = Proteome_Total_iBAQ_Norm > 0,
    .groups = "drop"
  )

# Abundance during proteomics dates = any non-zero metagenomic abundance on those dates
abundance_by_species <- species_abundance %>%
  group_by(Species_ID) %>%
  summarise(
    Species_Abundance_in_Proteomics_Dates = sum(Species_Abundance, na.rm = TRUE),
    Abundance_In_Proteomics_Dates         = Species_Abundance_in_Proteomics_Dates > 0,
    .groups = "drop"
  )

# Merge QC into species_meta; fill missing values safely
species_meta <- species_meta %>%
  left_join(proteome_by_species,  by = "Species_ID") %>%
  left_join(abundance_by_species, by = "Species_ID") %>%
  mutate(
    Proteome_Total_iBAQ_Norm              = ifelse(is.na(Proteome_Total_iBAQ_Norm), 0, Proteome_Total_iBAQ_Norm),
    Proteome_Recovered                    = ifelse(is.na(Proteome_Recovered), FALSE, Proteome_Recovered),
    Species_Abundance_in_Proteomics_Dates = ifelse(is.na(Species_Abundance_in_Proteomics_Dates), 0, Species_Abundance_in_Proteomics_Dates),
    Abundance_In_Proteomics_Dates         = ifelse(is.na(Abundance_In_Proteomics_Dates), FALSE, Abundance_In_Proteomics_Dates)
  )

############################################################
# F) Single plotting group: trophic group OR "Filtered"
############################################################

# Group_For_Plots is the ONE grouping variable used for figures/summaries.
# Species get their trophic group ONLY if ALL requirements are met:
#   (1) trophic info available (proxy for >=80% completeness)
#   (2) proteome recovered
#   (3) abundance present on proteomics sampling dates
# Otherwise: "Filtered"
species_meta <- species_meta %>%
  mutate(
    Group_For_Plots = ifelse(
      (TrophicInfo_Available %in% TRUE) &
        (Proteome_Recovered %in% TRUE) &
        (Abundance_In_Proteomics_Dates %in% TRUE),
      TrophicGroup,
      "Filtered"
    ),
    Group_For_Plots = factor(Group_For_Plots, levels = c("Oligo_like", "Intermediate", "Copio_like", "Filtered"))
  )

############################################################
# G) Propagate Group_For_Plots (+ relevant QC fields) to abundance + protein tables
#    WITHOUT creating duplicated columns
############################################################

# Define exactly which columns should be added downstream.
# We will only join columns that are NOT already present in the target table.
meta_cols_for_downstream <- species_meta %>%
  select(
    Species_ID,
    TrophicIndex, TrophicGroup, TrophicInfo_Available,
    Proteome_Recovered, Abundance_In_Proteomics_Dates,
    Group_For_Plots
  )

# --- G1) Abundance table ---
# species_abundance currently has only: Species_ID, Date, Species_Abundance
# so we can join all meta_cols_for_downstream safely.
species_abundance_annotated <- species_abundance %>%
  left_join(meta_cols_for_downstream, by = "Species_ID")

# --- G2) Protein table ---
# protein_table_final currently does not include these meta columns, so safe.
protein_table_final <- protein_table_final %>%
  left_join(meta_cols_for_downstream, by = "Species_ID")

############################################################
# H) Gene table + KEGG join + Group_For_Plots/QC (ALL species retained)
#    WITHOUT creating duplicated columns
############################################################

# Build gene table core (no species_meta columns yet)
gene_table <- gene_dat %>%
  select(
    Species_ID, Contig_ID, Protein_ID,
    KEGG_ID, KEGG_Annotation,
    COG_ID,  COG_Annotation,
    TIGR_ID, TIGR_Annotation,
    PFAM_ID, PFAM_Annotation,
    PHOBIUS_Annotation,
    Protein_Length, Protein_Sequence
  ) %>%
  distinct() %>%
  left_join(
    kegg_dat %>%
      select(
        KEGG_ID,
        KEGG_Category, KEGG_Subcategory_ID, KEGG_Subcategory,
        Transporter_Family, Substracte_Group, Substrate_Category
      ) %>%
      distinct(),
    by = "KEGG_ID"
  )

# Join species-level meta columns that are NOT already present in gene_table
# (defensive: avoids .x/.y if gene_dat ever changes in the future).
gene_table <- gene_table %>%
  left_join(
    meta_cols_for_downstream %>%
      select(Species_ID, TrophicIndex, TrophicGroup, TrophicInfo_Available,
             Proteome_Recovered, Abundance_In_Proteomics_Dates,
             Group_For_Plots),
    by = "Species_ID"
  )

############################################################
# I) Export final analysis tables
############################################################

write_csv(species_meta,              "Species_Meta_Table.csv")
write_csv(species_abundance_annotated, "Species_Abundance_Table.csv")
write_csv(protein_table_final,       "Species_Protein_Table.csv")
write_csv(gene_table,                "Species_Gene_Table.csv")

cat(
  "\nExport completed:\n",
  " - Species_Meta_Table.csv\n",
  " - Species_Abundance_Table.csv\n",
  " - Species_Protein_Table.csv\n",
  " - Species_Gene_Table.csv\n"
)

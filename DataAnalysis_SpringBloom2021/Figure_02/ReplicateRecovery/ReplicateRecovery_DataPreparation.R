################################################################################
#                                                                              #
#                             Replicate Recovery                               #
#                                                                              #
#   Author  : Cyrill Hofer                                                     #
#   Date    : 09.04.2025                                                       #
#   Purpose : Data Prep (median-centered by Date) + Species_Abundance in output#
#                                                                              #
################################################################################

rm(list = ls())

library(readr)
library(dplyr)

setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Figure_02/ReplicateRecovery/")

# ---------- Helper: median-centering by a grouping variable (Date) ----------
median_center_by_group <- function(df, group_col, value_col, out_col) {
  gcol <- rlang::ensym(group_col)
  vcol <- rlang::ensym(value_col)
  ocol <- rlang::ensym(out_col)
  
  per_group_med <- df %>%
    filter(!!vcol > 0) %>%
    group_by(!!gcol) %>%
    summarise(med_val = median(!!vcol, na.rm = TRUE), .groups = "drop")
  
  global_med <- median(per_group_med$med_val, na.rm = TRUE)
  
  df %>%
    left_join(per_group_med, by = rlang::as_name(gcol)) %>%
    mutate(!!ocol := (!!vcol) / med_val * global_med) %>%
    select(-med_val)
}

# ------------------------------- Load data -----------------------------------
Peptide_Dat <- read_csv(
  "../../../Datasets/72Samples_CostumDB/Integrated_Peptide_Table.csv",
  col_select = c(Date, Peptide_ID, Species_ID, Protein_ID, Peptide_Abundance,
                 Peptide_Presence, Phylum)
)

Protein_Dat <- read_csv(
  "../../../Datasets/72Samples_CostumDB/Integrated_Protein_Table.csv",
  col_select = c(Date, Protein_ID, Species_ID, Nr_Peptides, Nr_Tryptic_Peptides,
                 Protein_Abundance_IBAQ, Protein_Presence, Phylum, Class,
                 Order, Family, Genus, Species, Species_Abundance)
)

# remove contaminants
Peptide_Dat <- Peptide_Dat %>%
  filter(Phylum != "Y-FGCZCont")

Protein_Dat <- Protein_Dat %>%
  filter(Phylum != "Y-FGCZCont")

# ---------------------- Normalize (median-center by Date) --------------------
Peptide_Dat <- median_center_by_group(
  Peptide_Dat,
  group_col = Date,
  value_col = Peptide_Abundance,
  out_col   = Peptide_Abundance_Norm
)

Protein_Dat <- median_center_by_group(
  Protein_Dat,
  group_col = Date,
  value_col = Protein_Abundance_IBAQ,
  out_col   = Protein_Abundance_IBAQ_Norm
)

# ---- Build a per-Date x Species_ID lookup for Species_Abundance (robust) ----
species_abund_lookup <- Protein_Dat %>%
  group_by(Date, Species_ID) %>%
  summarise(
    Species_Abundance = if (all(is.na(Species_Abundance)))
      NA_real_ else median(Species_Abundance, na.rm = TRUE),
    .groups = "drop"
  )

# ---------------------- Long tables (only normalized values) -----------------
peptides_long <- Peptide_Dat %>%
  transmute(
    Date,
    ID        = Protein_ID,
    Level     = "Peptide",
    Abundance = Peptide_Abundance_Norm,    # normalized
    Presence  = Peptide_Presence,
    Species_ID,
    Phylum
  ) %>%
  # attach Species_Abundance by Date + Species_ID
  left_join(species_abund_lookup, by = c("Date", "Species_ID"))

proteins_long <- Protein_Dat %>%
  transmute(
    Date,
    ID        = Protein_ID,
    Level     = "Protein",
    Abundance = Protein_Abundance_IBAQ_Norm,  # normalized
    Presence  = Protein_Presence,
    Species_ID,
    Phylum,
    Species_Abundance                           # pass through existing
  )

Combined_Data <- bind_rows(peptides_long, proteins_long)

# --------------------------------- Export ------------------------------------
write_csv(Combined_Data, "Proteomic_Data_Overview_Level_Long.csv")

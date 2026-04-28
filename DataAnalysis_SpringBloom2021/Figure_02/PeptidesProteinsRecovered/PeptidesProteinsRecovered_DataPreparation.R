################################################################################
# Proteome specificity + DIA-NN resolution annotation
# + identity-aware peptide specificity (UNFOLDING Protein.Group ONLY)
# Cyrill Hofer | 13.01.2026
#
# Reviewer-safe definition used here:
#   - "Protein groups" are taken exactly as reported in DIA-NN `Protein.Group`.
#   - If `Protein.Group` contains multiple accessions separated by ';',
#     we treat these as the alternative members of that protein group and
#     compute taxonomy resolution from this unfolded list.
#   - We DO NOT use `Protein.Ids` to expand membership (by design), to stay
#     consistent with DIA-NN’s reported group representation.
#
# Outputs (written to working directory):
#  1) diann_annotated_withProteinAndPeptideResolution.csv
#     - full DIA-NN table with:
#         ProteinGroup_Resolution, Species_Specific_ProteinGroup
#         Peptide_Resolution_IdentityAware, Species_Specific_Peptide
#         n_ProteinGroups, n_Proteins_Unfolded (per peptide)
#
#  2) ProteinGroupResolution_LUT.csv
#     - one row per DIA-NN Protein.Group
#
#  3) PeptideResolution_IdentityAware_LUT.csv
#     - one row per peptide (Stripped.Sequence)
#
#  4) PeptideProteinDistribution_withResolution.csv
#     - stage summary (Raw / Q<1% / All<1% / 2+ pep) including species-specific counts
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# ------------------------- USER SETTINGS (edit paths) --------------------------
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Figure_02/PeptidesProteinsRecovered/")

PATH_DIANN   <- "../../../Data/Proteomics_Data/Raw_Data/72Samples_p34548_SpringBloom2021_95MagSpeciesDB_15JUL24_20240726_Database_Run07JAN25/diann-report.tsv"
PATH_SPECIES <- "../../../Data/Species_Data/Final_Data/Species_Table.csv"

# ------------------------------ 1) READ DATA ---------------------------------
diann <- read_tsv(PATH_DIANN, show_col_types = FALSE)

Species_Dat <- read_csv(
  PATH_SPECIES,
  col_select = c(Species_ID, Phylum, Class, Order, Family, Genus, Species),
  show_col_types = FALSE
) %>% distinct()

# ---------------------- Remove contaminants (Y-FGCZCont) ----------------------
# NOTE: filtering contaminants is orthogonal to resolution logic; we keep it here
# to prevent contaminant identifiers from inflating ambiguity.

diann <- diann %>%
  mutate(
    is_contam_YFGCZ = str_detect(`Protein.Group`, "Y-FGCZCont") |
      str_detect(`Protein.Ids`,   "Y-FGCZCont")
  )

# optional: remove DIA-NN decoys if column exists
if ("Decoy" %in% colnames(diann)) {
  diann <- diann %>%
    mutate(is_decoy = as.logical(Decoy))
} else {
  diann <- diann %>% mutate(is_decoy = FALSE)
}

# report + filter
diann %>%
  summarise(
    n_total        = n(),
    n_YFGCZcont    = sum(is_contam_YFGCZ, na.rm = TRUE),
    n_decoy        = sum(is_decoy, na.rm = TRUE),
    n_removed_any  = sum(is_contam_YFGCZ | is_decoy, na.rm = TRUE),
    n_remaining    = sum(!(is_contam_YFGCZ | is_decoy), na.rm = TRUE)
  ) %>% print()

diann <- diann %>%
  filter(!(is_contam_YFGCZ | is_decoy)) %>%
  select(-is_contam_YFGCZ, -is_decoy)

# ==============================================================================
# 2) PROTEIN-GROUP RESOLUTION (UNFOLD Protein.Group ONLY)
# ------------------------------------------------------------------------------
# Goal:
#   For each DIA-NN Protein.Group, determine the finest taxonomic rank shared
#   by all unfolded members of Protein.Group (LCA-style resolution).
#
# Implementation:
#   - Build unique list of Protein.Group values.
#   - Unfold Protein.Group by splitting on ';' into GroupMember.
#   - Map GroupMember -> Species_ID (based on MAG naming convention).
#   - Join taxonomy.
#   - Classify resolution using n_distinct at each rank.
#
# Conservative choice:
#   - If ANY unfolded member lacks taxonomy (Phylum is NA), assign "No_Match".
#     (This is intentionally conservative and will reduce "species-specific"
#      calls when any member cannot be mapped.)
# ==============================================================================

protein_groups <- diann %>%
  distinct(Protein.Group) %>%
  filter(!is.na(Protein.Group), Protein.Group != "")

group_members_pg <- protein_groups %>%
  # Unfold Protein.Group itself (semicolon-separated members)
  mutate(GroupMember = Protein.Group) %>%
  separate_rows(GroupMember, sep = ";") %>%
  mutate(
    GroupMember = str_trim(GroupMember),
    # Extract Species_ID only for MAG-style accessions; else NA
    Species_ID = if_else(
      str_detect(GroupMember, "\\.MAG_"),
      str_extract(GroupMember, "^[^-]+"),
      NA_character_
    )
  ) %>%
  filter(!is.na(GroupMember), GroupMember != "") %>%
  left_join(Species_Dat, by = "Species_ID")

protein_group_resolution_lut <- group_members_pg %>%
  group_by(Protein.Group) %>%
  summarise(
    ProteinGroup_Resolution = case_when(
      any(is.na(Phylum)) ~ "No_Match",
      n_distinct(Species, na.rm = TRUE) == 1 ~ "Species_Resolution",
      n_distinct(Genus,   na.rm = TRUE) == 1 ~ "Genus_Resolution",
      n_distinct(Family,  na.rm = TRUE) == 1 ~ "Family_Resolution",
      n_distinct(Order,   na.rm = TRUE) == 1 ~ "Order_Resolution",
      n_distinct(Class,   na.rm = TRUE) == 1 ~ "Class_Resolution",
      n_distinct(Phylum,  na.rm = TRUE) == 1 ~ "Phylum_Resolution",
      TRUE ~ "No_Match"
    ),
    Species_Specific_ProteinGroup = (ProteinGroup_Resolution == "Species_Resolution"),
    .groups = "drop"
  )

# Join protein-group resolution onto full DIA-NN table
diann_pg <- diann %>%
  left_join(protein_group_resolution_lut, by = "Protein.Group")

# ==============================================================================
# 3) IDENTITY-AWARE PEPTIDE SPECIFICITY (UNFOLD Protein.Group ONLY)
# ------------------------------------------------------------------------------
# Goal:
#   For each peptide (Stripped.Sequence), compute identity-aware specificity by:
#     - collecting all Protein.Groups the peptide maps to,
#     - unfolding each Protein.Group (semicolon-separated members),
#     - evaluating taxonomy across ALL possible members (LCA-style).
#
# Note:
#   This is conservative and reflects identification ambiguity, not expression.
# ==============================================================================

# 3a) peptide -> protein-group mapping (unique)
pep_pg <- diann %>%
  select(Stripped.Sequence, Protein.Group) %>%
  filter(!is.na(Stripped.Sequence), !is.na(Protein.Group),
         Stripped.Sequence != "", Protein.Group != "") %>%
  distinct()

# 3b) Expand each peptide to all unfolded Protein.Group members
pep_member_tax <- pep_pg %>%
  mutate(GroupMember = Protein.Group) %>%
  separate_rows(GroupMember, sep = ";") %>%
  mutate(
    GroupMember = str_trim(GroupMember),
    Species_ID  = if_else(
      str_detect(GroupMember, "\\.MAG_"),
      str_extract(GroupMember, "^[^-]+"),
      NA_character_
    )
  ) %>%
  distinct(Stripped.Sequence, Protein.Group, GroupMember, Species_ID) %>%
  left_join(Species_Dat, by = "Species_ID")

# 3c) Summarise peptide LCA-style resolution
peptide_specificity_tbl <- pep_member_tax %>%
  group_by(Stripped.Sequence) %>%
  summarise(
    # how many protein groups this peptide maps to
    n_ProteinGroups = n_distinct(Protein.Group),
    # how many unfolded members across those groups
    n_Proteins_Unfolded = n_distinct(GroupMember),
    
    # uniqueness at each rank (identity-aware)
    n_Species = n_distinct(Species, na.rm = TRUE),
    n_Genus   = n_distinct(Genus,   na.rm = TRUE),
    n_Family  = n_distinct(Family,  na.rm = TRUE),
    n_Order   = n_distinct(Order,   na.rm = TRUE),
    n_Class   = n_distinct(Class,   na.rm = TRUE),
    n_Phylum  = n_distinct(Phylum,  na.rm = TRUE),
    
    # conservative: if any member lacks taxonomy -> No_Match
    any_missing_taxonomy = any(is.na(Phylum)),
    
    Peptide_Resolution_IdentityAware = case_when(
      any_missing_taxonomy ~ "No_Match",
      n_Species == 1 ~ "Species_Resolution",
      n_Genus   == 1 ~ "Genus_Resolution",
      n_Family  == 1 ~ "Family_Resolution",
      n_Order   == 1 ~ "Order_Resolution",
      n_Class   == 1 ~ "Class_Resolution",
      n_Phylum  == 1 ~ "Phylum_Resolution",
      TRUE           ~ "No_Match"
    ),
    
    Species_Specific_Peptide = (Peptide_Resolution_IdentityAware == "Species_Resolution"),
    .groups = "drop"
  )

# Add peptide specificity to DIA-NN
diann_full <- diann_pg %>%
  left_join(peptide_specificity_tbl, by = "Stripped.Sequence")

# ==============================================================================
# 4) FILTERING STAGES (unchanged)
# ==============================================================================

f1 <- diann_full %>% filter(Q.Value < 0.01)

f2 <- diann_full %>% filter(Q.Value < 0.01,
                            Lib.Q.Value < 0.01,
                            PG.Q.Value < 0.01)

f3 <- diann_full %>%
  filter(Q.Value < 0.01,
         Lib.Q.Value < 0.01,
         PG.Q.Value < 0.01) %>%
  group_by(Protein.Group) %>%
  filter(n_distinct(Stripped.Sequence) >= 2) %>%
  ungroup()

stages <- list(Raw = diann_full, `Q<1%` = f1, `All<1%` = f2, `2+ pep` = f3)

# ==============================================================================
# 5) STAGE SUMMARY TABLE (unchanged)
# ==============================================================================

summary_tbl <- purrr::imap_dfr(stages, ~{
  df <- .x
  pep_ids <- df %>% distinct(Stripped.Sequence, Proteotypic)
  
  tibble(
    Stage = .y,
    
    # totals
    Peptides = n_distinct(df$Stripped.Sequence),
    Proteins = n_distinct(df$Protein.Group),
    
    # DIA-NN proteotypic evidence
    ProteotypicPeptides = sum(pep_ids$Proteotypic == 1, na.rm = TRUE),
    ProteotypicProteins = df %>% filter(Proteotypic == 1) %>% distinct(Protein.Group) %>% nrow(),
    
    # protein-group species-specific
    SpeciesSpecificProteinGroups = df %>%
      distinct(Protein.Group, Species_Specific_ProteinGroup) %>%
      filter(Species_Specific_ProteinGroup %in% TRUE) %>%
      nrow(),
    
    # peptide species-specific (identity-aware)
    SpeciesSpecificPeptides = df %>%
      distinct(Stripped.Sequence, Species_Specific_Peptide) %>%
      filter(Species_Specific_Peptide %in% TRUE) %>%
      nrow()
  )
}) %>%
  mutate(Stage = factor(Stage, levels = c("Raw", "Q<1%", "All<1%", "2+ pep")))

# ==============================================================================
# 6) SUMMARY PRINT
# ==============================================================================

cat("\n==================== SUMMARY BEFORE EXPORT ====================\n")
cat("DIA-NN rows (filtered):             ", nrow(diann), "\n")
cat("Unique Protein.Groups (DIA-NN):     ", n_distinct(diann$Protein.Group), "\n")
cat("Unique peptides (DIA-NN):           ", n_distinct(diann$Stripped.Sequence), "\n")
cat("Protein.Groups with resolution:     ", nrow(protein_group_resolution_lut), "\n")
cat("Peptides with identity-aware resolution: ", nrow(peptide_specificity_tbl), "\n\n")
print(summary_tbl)
cat("===============================================================\n\n")

# ==============================================================================
# 7) EXPORT OUTPUTS
# ==============================================================================

write_csv(diann_full, "diann_annotated_withProteinAndPeptideResolution.csv")
write_csv(protein_group_resolution_lut, "ProteinGroupResolution_LUT.csv")
write_csv(peptide_specificity_tbl, "PeptideResolution_IdentityAware_LUT.csv")
write_csv(summary_tbl, "PeptideProteinDistribution_withResolution.csv")

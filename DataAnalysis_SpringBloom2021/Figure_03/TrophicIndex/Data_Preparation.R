# =============================================================================
# Trophic Index Data Preparation
#
# Build:
#   1) Species_GenomicProperties.csv
#   2) Species_GenomicProperties_Completed.csv
#   3) Species_GenomicProperties_ObservedVsEstimated.csv   <-- NEW (combined)
#
# Combined file contains:
#   - Completeness + Contamination (+ Estimated_Genome_Size, Total_CDS)
#   - For each gene-count-like feature: observed + estimated (side-by-side)
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
})

# -----------------------------
# 0) Paths
# -----------------------------
WORKDIR      <- "/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Figure_03/TrophicIndex/"
PATH_SPECIES <- "../../../Data/Species_Data/Final_Data/Species_Table.csv"
PATH_MAG     <- "../../../Data/MAG_Data/Final_Data/MAG_Table.csv"

OUTFILE_OBS   <- "Species_GenomicProperties.csv"
OUTFILE_COMP  <- "Species_GenomicProperties_Completed.csv"
OUTFILE_BOTH  <- "Species_GenomicProperties_ObservedVsEstimated.csv"

setwd(WORKDIR)

# -----------------------------
# 1) Read inputs
# -----------------------------
species_table <- read_csv(PATH_SPECIES, show_col_types = FALSE)
mag_table     <- read_csv(PATH_MAG, show_col_types = FALSE)

# -----------------------------
# 2) Keep only representative, high-quality genomes
# -----------------------------
mag_keep <- mag_table %>%
  filter(Species_ID == MAG_ID, Completeness >= 80) %>%
  distinct(Species_ID, .keep_all = TRUE)

keep_species <- mag_keep$Species_ID

# -----------------------------
# 3) Taxonomy (one row per Species_ID)
# -----------------------------
tax_keep <- species_table %>%
  filter(Species_ID %in% keep_species) %>%
  select(Species_ID, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  distinct(Species_ID, .keep_all = TRUE)

# -----------------------------
# 4) Accession parsing
# -----------------------------
acc_keep <- mag_keep %>%
  transmute(
    Species_ID,
    Accession = str_extract(Accession_Number, "^GCA_[0-9]+")
  )

# -----------------------------
# 5) Observed features
# -----------------------------
feat_keep <- mag_keep %>%
  transmute(
    Species_ID,
    GC = GC_Content,
    GENOME_LENGTH = Nr_Bases,
    CODING_DENSITY = Coding_Density,
    MEAN_INTERGENIC_SPACER = Mean_Intergenic_Spacer,
    TOTAL_TRNA = Total_tRNA,
    NUM_5S  = Nr_5s,
    NUM_16S = Nr_16s,
    NUM_23S = Nr_23s,
    NUM_SIGMA = Nr_Sigma,
    NUM_SIGNALTRANS = Nr_Signaltrans
  ) %>%
  mutate(
    MAX_RRNA = pmax(NUM_5S, NUM_16S, NUM_23S, na.rm = TRUE)
  )

# -----------------------------
# 6) MAG meta (include completeness/contamination for ALL downstream outputs)
# -----------------------------
mag_meta <- mag_table %>%
  filter(Species_ID == MAG_ID) %>%
  select(
    Species_ID,
    Completeness,
    Contamination,
    Estimated_Genome_Size,
    Total_CDS
  ) %>%
  distinct(Species_ID, .keep_all = TRUE)

# -----------------------------
# 7) final_dat (observed)  [now includes completeness/contamination]
# -----------------------------
final_dat <- tax_keep %>%
  left_join(acc_keep, by = "Species_ID") %>%
  left_join(mag_meta, by = "Species_ID") %>%
  left_join(feat_keep, by = "Species_ID") %>%
  select(
    Species_ID,
    Domain, Phylum, Class, Order, Family, Genus, Species,
    Accession,
    Completeness, Contamination, Estimated_Genome_Size, Total_CDS,
    GC, GENOME_LENGTH,
    TOTAL_TRNA, NUM_5S, NUM_16S, NUM_23S, MAX_RRNA,
    CODING_DENSITY, MEAN_INTERGENIC_SPACER, NUM_SIGMA, NUM_SIGNALTRANS
  ) %>%
  filter(
    Species != "Planktothrix agardhii",
    Order   != "Rickettsiales",
    Phylum  != "Patescibacteria",
    Domain  != "Archaea"
  )

stopifnot(!anyDuplicated(final_dat$Species_ID))

write_csv(final_dat, OUTFILE_OBS)
cat("\nWrote:", OUTFILE_OBS, "\n")
glimpse(final_dat)

# -----------------------------
# 8) final_dat_completed (completion-corrected; keeps observed columns too)
# -----------------------------
final_dat_completed <- final_dat %>%
  mutate(
    Estimated_Genome_Size = as.numeric(Estimated_Genome_Size),
    Total_CDS             = as.numeric(Total_CDS),
    Completeness          = as.numeric(Completeness),
    Contamination         = as.numeric(Contamination),
    
    # missing genome length (bp)
    delta_length = pmax(0, Estimated_Genome_Size - GENOME_LENGTH),
    
    # gene density and missing CDS
    gene_density = Total_CDS / GENOME_LENGTH,
    Missing_CDS  = gene_density * delta_length,
    
    # fractions among recovered genes
    frac_trna        = ifelse(Total_CDS > 0, TOTAL_TRNA      / Total_CDS, NA_real_),
    frac_sigma       = ifelse(Total_CDS > 0, NUM_SIGMA       / Total_CDS, NA_real_),
    frac_signaltrans = ifelse(Total_CDS > 0, NUM_SIGNALTRANS / Total_CDS, NA_real_),
    
    frac_5S          = ifelse(Total_CDS > 0, NUM_5S          / Total_CDS, NA_real_),
    frac_16S         = ifelse(Total_CDS > 0, NUM_16S         / Total_CDS, NA_real_),
    frac_23S         = ifelse(Total_CDS > 0, NUM_23S         / Total_CDS, NA_real_),
    
    # expected missing features
    Missing_TRNA        = Missing_CDS * frac_trna,
    Missing_SIGMA       = Missing_CDS * frac_sigma,
    Missing_SIGNALTRANS = Missing_CDS * frac_signaltrans,
    
    Missing_5S          = Missing_CDS * frac_5S,
    Missing_16S         = Missing_CDS * frac_16S,
    Missing_23S         = Missing_CDS * frac_23S,
    
    # completion-corrected totals
    Est_TOTAL_TRNA      = TOTAL_TRNA      + Missing_TRNA,
    Est_NUM_SIGMA       = NUM_SIGMA       + Missing_SIGMA,
    Est_NUM_SIGNALTRANS = NUM_SIGNALTRANS + Missing_SIGNALTRANS,
    
    Est_NUM_5S          = NUM_5S          + Missing_5S,
    Est_NUM_16S         = NUM_16S         + Missing_16S,
    Est_NUM_23S         = NUM_23S         + Missing_23S,
    Est_MAX_RRNA        = pmax(Est_NUM_5S, Est_NUM_16S, Est_NUM_23S, na.rm = TRUE)
  )

stopifnot(!anyDuplicated(final_dat_completed$Species_ID))

write_csv(final_dat_completed, OUTFILE_COMP)
cat("\nWrote:", OUTFILE_COMP, "\n")
glimpse(final_dat_completed)

# -----------------------------
# 9) Combined "Observed vs Estimated" dataset (side-by-side)
#     - Keeps Completeness/Contamination
#     - Has both observed + estimated columns for corrected features
# -----------------------------
final_dat_both <- final_dat %>%
  left_join(
    final_dat_completed %>%
      select(
        Species_ID,
        Est_TOTAL_TRNA, Est_NUM_5S, Est_NUM_16S, Est_NUM_23S, Est_MAX_RRNA,
        Est_NUM_SIGMA, Est_NUM_SIGNALTRANS,
        delta_length, gene_density, Missing_CDS
      ),
    by = "Species_ID"
  ) %>%
  # make it explicit these are "observed"
  rename(
    Obs_TOTAL_TRNA      = TOTAL_TRNA,
    Obs_NUM_5S          = NUM_5S,
    Obs_NUM_16S         = NUM_16S,
    Obs_NUM_23S         = NUM_23S,
    Obs_MAX_RRNA        = MAX_RRNA,
    Obs_NUM_SIGMA       = NUM_SIGMA,
    Obs_NUM_SIGNALTRANS = NUM_SIGNALTRANS
  ) %>%
  select(
    Species_ID,
    Domain, Phylum, Class, Order, Family, Genus, Species,
    Accession,
    Completeness, Contamination, Estimated_Genome_Size, Total_CDS,
    GC, GENOME_LENGTH, CODING_DENSITY, MEAN_INTERGENIC_SPACER,
    
    # side-by-side observed vs estimated
    Obs_TOTAL_TRNA,      Est_TOTAL_TRNA,
    Obs_NUM_5S,          Est_NUM_5S,
    Obs_NUM_16S,         Est_NUM_16S,
    Obs_NUM_23S,         Est_NUM_23S,
    Obs_MAX_RRNA,        Est_MAX_RRNA,
    Obs_NUM_SIGMA,       Est_NUM_SIGMA,
    Obs_NUM_SIGNALTRANS, Est_NUM_SIGNALTRANS,
    
    # diagnostics (helpful in review / methods)
    delta_length, gene_density, Missing_CDS
  )

stopifnot(!anyDuplicated(final_dat_both$Species_ID))

write_csv(final_dat_both, OUTFILE_BOTH)
cat("\nWrote:", OUTFILE_BOTH, "\n")
glimpse(final_dat_both)
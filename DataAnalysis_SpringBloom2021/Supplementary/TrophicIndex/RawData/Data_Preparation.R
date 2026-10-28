# =============================================================================
# Indicator species — Genomic Properties ONLY (no COG profiles)
# Output:
#   ../IndicatorSpecies_GenomicProperties.csv
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
})

# ---- Working directory (raw indicator inputs live here) ----
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Figure_03/TrophicIndex/RawData/")

# ---- Read inputs ----
species        <- read.delim("Genome_Info.tsv", header = FALSE)          # species metadata
features       <- read.delim("GenomicFeatures.tsv")                      # genome-level features
trophic_status <- read_excel("ReferenceList_TrophicStrategies.xlsx")     # Accession -> Classification

# ---- Standardize species table + accession join key ----
# Species_ID examples (as in your notes):
#   - Genome_CP000655.1  -> Accession = CP000655.1
#   - GCA_965194435.1_*  -> Accession = GCA_965194435  (drop version/suffix)
species <- species %>%
  rename(
    Species_ID    = V1,
    Contig_ID     = V2,
    Taxonomy      = V3,
    Genome_Status = V4
  ) %>%
  mutate(
    Accession = case_when(
      str_detect(Species_ID, "^Genome_") ~ str_remove(Species_ID, "^Genome_"),
      str_detect(Species_ID, "^GCA_")    ~ str_extract(Species_ID, "^GCA_[0-9]+"),
      TRUE                               ~ NA_character_
    )
  )

# ---- Trophic strategy lookup (keep only required columns) ----
trophic_status <- trophic_status %>%
  select(Accession, Classification) %>%
  distinct()

# ---- Standardize features table join key ----
# (your previous pipeline used: rename(Species_ID = FILENAME))
features <- features %>%
  rename(Species_ID = FILENAME)

# ---- Select/rename genomic features to match your downstream workflow ----
# IMPORTANT: keep these names stable (they are used later in PCA scripts)
features_selected <- features %>%
  transmute(
    Species_ID,
    GC               = GC,             # keep as-is if already GC
    GENOME_LENGTH    = GENOME_LENGTH,  # keep as-is if already GENOME_LENGTH
    TOTAL_TRNA       = TOTAL_TRNA,
    NUM_5S           = NUM_5S,
    NUM_16S          = NUM_16S,
    NUM_23S          = NUM_23S,
    CODING_DENSITY   = CODING_DENSITY,
    MEAN_INTERGENIC_SPACER = MEAN_INTERGENIC_SPACER,
    NUM_SIGMA        = NUM_SIGMA,
    NUM_SIGNALTRANS  = NUM_SIGNALTRANS,
    
    # Derived feature: maximum rRNA copy number across loci
    MAX_RRNA = pmax(NUM_5S, NUM_16S, NUM_23S, na.rm = TRUE)
  ) %>%
  distinct()

# ---- Build final indicator training table (1 row per Species_ID) ----
final_dat <- species %>%
  select(Species_ID, Taxonomy, Accession) %>%
  distinct() %>%
  left_join(trophic_status, by = "Accession") %>%
  left_join(features_selected, by = "Species_ID")

# ---- Minimal sanity checks (fail loudly if something is off) ----
cat("\nRows (species):", nrow(final_dat), "\n")
cat("Missing Classification:", sum(is.na(final_dat$Classification)), "\n")
cat("Missing Accession:", sum(is.na(final_dat$Accession)), "\n")

# ---- Export ----
write_csv(final_dat, "../IndicatorSpecies_GenomicProperties.csv")
cat("\nWrote: ../IndicatorSpecies_GenomicProperties.csv\n")

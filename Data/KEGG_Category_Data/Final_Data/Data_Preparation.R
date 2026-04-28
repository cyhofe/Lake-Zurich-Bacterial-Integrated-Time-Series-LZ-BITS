# ============================================================
# KEGG BRITE – Global + Transporter family tables (MINIMAL)
# with expanded BRITE path levels for transporters
# ============================================================

rm(list = ls())

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/KEGG_Category_Data/Final_Data/")

infile <- "../Raw_Data/kegg_brite_master_long.tsv"

out_global <- "KEGG_Functional_Categories.csv"
out_trans  <- "KEGG_Transporter_Categories.csv"

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
raw_dat <- read_delim(
  infile,
  delim = "\t",
  col_types = cols()
)

# ------------------------------------------------------------
# 1) GLOBAL functional categories
# ------------------------------------------------------------
KEGG_Functional_Categories <- raw_dat %>%
  select(KO, major_group, brite_id, brite_name) %>%
  distinct() %>%
  arrange(major_group, brite_name, KO)

write_csv(KEGG_Functional_Categories, out_global)

# ------------------------------------------------------------
# 2) TRANSPORTERS ONLY – clean hierarchy extraction
# ------------------------------------------------------------

KEGG_Transporter_Categories <- raw_dat %>%
  filter(brite_id == "ko02000") %>%
  mutate(
    # split path into list-column
    path_levels = str_split(path, " > "),
    
    # transporter family (mechanism)
    Transporter_Family = case_when(
      str_detect(tolower(path), "phosphotransferase system|\\bpts\\b") ~ "PTS",
      str_detect(tolower(path), "major facilitator superfamily|\\bmfs\\b") ~ "MFS",
      str_detect(tolower(path), "\\bslc\\b|solute carrier") ~ "SLC",
      str_detect(tolower(path), "abc transporters|atp-binding cassette") ~ "ABC",
      TRUE ~ "Other"
    ),
    
    # extract meaningful hierarchy levels safely
    Path_L1 = map_chr(path_levels, ~ .x[1] %||% NA_character_),
    Path_L2 = map_chr(path_levels, ~ .x[2] %||% NA_character_),
    Path_L3 = map_chr(path_levels, ~ .x[3] %||% NA_character_),
    Path_L4 = map_chr(path_levels, ~ .x[4] %||% NA_character_),
    Path_L5 = map_chr(path_levels, ~ .x[5] %||% NA_character_)
  ) %>%
  select(
    KO,
    Transporter_Family,
    major_group,
    brite_id,
    brite_name,
    Path_L1,
    Path_L2,
    Path_L3,
    Path_L4,
    Path_L5,
    path
  ) %>%
  distinct() %>%
  arrange(Transporter_Family, KO)

write_csv(KEGG_Transporter_Categories, out_trans)

# ------------------------------------------------------------
# QC prints (short)
# ------------------------------------------------------------
cat("\nWritten files:\n")
cat(" -", out_global, "\n")
cat(" -", out_trans, "\n")

cat("\nTransporter family counts (unique KOs):\n")
print(
  KEGG_Transporter_Categories %>%
    distinct(KO, Transporter_Family) %>%
    count(Transporter_Family, sort = TRUE)
)

cat("\nMax BRITE path depth in transporters:\n")
print(
  KEGG_Transporter_Categories %>%
    transmute(depth = rowSums(!is.na(select(., starts_with("path_L"))))) %>%
    summarise(max_depth = max(depth))
)

cat("\n✅ Done.\n")



unique(KEGG_Transporter_Categories$Path_L4)

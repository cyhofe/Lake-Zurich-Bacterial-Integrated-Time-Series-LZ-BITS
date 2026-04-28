# ============================================================
# KEGG categories – SIMPLE summary statistics
# ============================================================

rm(list = ls())

library(dplyr)
library(readr)

setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/KEGG_Category_Data/Final_Data/")

# ------------------------------------------------------------
# Read and prepare data
# ------------------------------------------------------------
kegg_categories <- read_csv("KEGG_Functional_Categories.csv", show_col_types = FALSE) %>%
  rename(
    KEGG_ID = KO,
    KEGG_Category = major_group,
    KEGG_Subcategory_ID = brite_id,
    KEGG_Subcategory = brite_name
  ) %>%
  distinct()

transporter_categories <- read_csv("KEGG_Transporter_Categories_Curated.csv", show_col_types = FALSE) %>%
  select(KO, Transporter_Family, Substracte_Group, Substrate_Category) %>%
  rename(KEGG_ID = KO) %>%
  distinct()

# ============================================================
# GLOBAL KEGG CATEGORIES
# ============================================================

# KEGG IDs per major category
global_ids_per_category <- kegg_categories %>%
  group_by(KEGG_Category) %>%
  summarise(n_KEGG_IDs = n_distinct(KEGG_ID), .groups = "drop")

# KEGG IDs in one vs multiple categories
global_kegg_membership <- kegg_categories %>%
  distinct(KEGG_ID, KEGG_Category) %>%
  count(KEGG_ID, name = "n_categories")

global_kegg_membership_summary <- global_kegg_membership %>%
  summarise(
    total_KEGG_IDs = n(),
    unique_to_one_category = sum(n_categories == 1),
    in_multiple_categories = sum(n_categories > 1)
  )

# ============================================================
# TRANSPORTER CATEGORIES
# ============================================================

# KEGG IDs per transporter family
trans_ids_per_family <- transporter_categories %>%
  group_by(Transporter_Family) %>%
  summarise(n_KEGG_IDs = n_distinct(KEGG_ID), .groups = "drop")

# KEGG IDs per substrate category
trans_ids_per_substrate_cat <- transporter_categories %>%
  group_by(Substrate_Category) %>%
  summarise(n_KEGG_IDs = n_distinct(KEGG_ID), .groups = "drop")

# KEGG IDs per substrate group
trans_ids_per_substrate_group <- transporter_categories %>%
  group_by(Substracte_Group) %>%
  summarise(n_KEGG_IDs = n_distinct(KEGG_ID), .groups = "drop")

# KEGG IDs in one vs multiple categories (per column)
trans_kegg_membership <- transporter_categories %>%
  distinct(KEGG_ID, Transporter_Family, Substracte_Group, Substrate_Category) %>%
  pivot_longer(
    cols = c(Transporter_Family, Substracte_Group, Substrate_Category),
    names_to = "category_type",
    values_to = "category"
  ) %>%
  filter(!is.na(category)) %>%
  distinct(KEGG_ID, category_type, category) %>%
  count(category_type, KEGG_ID, name = "n_categories")

trans_kegg_membership_summary <- trans_kegg_membership %>%
  group_by(category_type) %>%
  summarise(
    total_KEGG_IDs = n(),
    unique_to_one_category = sum(n_categories == 1),
    in_multiple_categories = sum(n_categories > 1),
    .groups = "drop"
  )

# ------------------------------------------------------------
# View results
# ------------------------------------------------------------
global_ids_per_category
global_kegg_membership_summary

trans_ids_per_family
trans_ids_per_substrate_cat
trans_ids_per_substrate_group
trans_kegg_membership_summary

# ------------------------------------------------------------
# Collapse global KEGG categories to one row per KEGG_ID
# ------------------------------------------------------------
kegg_collapsed <- kegg_categories %>%
  group_by(KEGG_ID) %>%
  summarise(
    KEGG_Category = paste(sort(unique(KEGG_Category)), collapse = "; "),
    KEGG_Subcategory_ID = paste(sort(unique(KEGG_Subcategory_ID)), collapse = "; "),
    KEGG_Subcategory = paste(sort(unique(KEGG_Subcategory)), collapse = "; "),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Collapse transporter categories to one row per KEGG_ID
# ------------------------------------------------------------
transporter_collapsed <- transporter_categories %>%
  group_by(KEGG_ID) %>%
  summarise(
    Transporter_Family = paste(sort(unique(Transporter_Family)), collapse = "; "),
    Substracte_Group = paste(sort(unique(Substracte_Group)), collapse = "; "),
    Substrate_Category = paste(sort(unique(Substrate_Category)), collapse = "; "),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Final join (ONE row per KEGG_ID)
# ------------------------------------------------------------
final_kegg_table <- full_join(
  kegg_collapsed,
  transporter_collapsed,
  by = "KEGG_ID"
)

# View
final_kegg_table

# KEGG IDs in each table
kegg_ids_global <- unique(kegg_categories$KEGG_ID)
kegg_ids_trans  <- unique(transporter_categories$KEGG_ID)
kegg_ids_final  <- unique(final_kegg_table$KEGG_ID)

# Check: all global KEGG IDs are present
all(kegg_ids_global %in% kegg_ids_final)

# Check: all transporter KEGG IDs are present
all(kegg_ids_trans %in% kegg_ids_final)

length(kegg_ids_global)
length(kegg_ids_trans)
length(kegg_ids_final)
length(kegg_ids_final) ==
  length(unique(c(kegg_ids_global, kegg_ids_trans)))

# ------------------------------------------------------------
# Export final combined KEGG table
# ------------------------------------------------------------
write.csv(
  final_kegg_table,
  file = "KEGG_Category_Table.csv",
  row.names = FALSE
)

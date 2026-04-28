# Clear environment
rm(list = ls())

library(dplyr)
library(readr)

# Load data
old <- read_csv("Gene_Table.csv")
new <- read_csv("Gene_Table_2.csv")

cat("===== BASIC ROW COMPARISON =====\n")
cat("Old row count: ", nrow(old), "\n")
cat("New row count: ", nrow(new), "\n\n")

cat("===== UNIQUE PROTEIN IDS =====\n")
cat("Old unique Protein_IDs: ", n_distinct(old$Protein_ID), "\n")
cat("New unique Protein_IDs: ", n_distinct(new$Protein_ID), "\n")

missing_in_new <- setdiff(old$Protein_ID, new$Protein_ID)
missing_in_old <- setdiff(new$Protein_ID, old$Protein_ID)

cat("Protein_IDs missing in new table: ", length(missing_in_new), "\n")
cat("Protein_IDs missing in old table: ", length(missing_in_old), "\n\n")

# KEGG comparison --------------------------------------------------------------

cat("===== KEGG ID COMPARISON =====\n")
cat("Old KEGG_ID distinct: ", n_distinct(old$KEGG_ID, na.rm = TRUE), "\n")
cat("New KEGG_ID distinct: ", n_distinct(new$KEGG_ID, na.rm = TRUE), "\n")

missing_kegg_new <- setdiff(na.omit(old$KEGG_ID), na.omit(new$KEGG_ID))
missing_kegg_old <- setdiff(na.omit(new$KEGG_ID), na.omit(old$KEGG_ID))

cat("KEGG IDs missing in new: ", length(missing_kegg_new), "\n")
cat("KEGG IDs missing in old: ", length(missing_kegg_old), "\n\n")

# COG comparison ---------------------------------------------------------------

cat("===== COG ID COMPARISON =====\n")
cat("Old COG_ID distinct: ", n_distinct(old$COG_ID, na.rm = TRUE), "\n")
cat("New COG_ID distinct: ", n_distinct(new$COG_ID, na.rm = TRUE), "\n")

missing_cog_new <- setdiff(na.omit(old$COG_ID), na.omit(new$COG_ID))
missing_cog_old <- setdiff(na.omit(new$COG_ID), na.omit(old$COG_ID))

cat("COG IDs missing in new: ", length(missing_cog_new), "\n")
cat("COG IDs missing in old: ", length(missing_cog_old), "\n\n")

# Check for duplicated Protein_IDs ---------------------------------------------

cat("===== DUPLICATION CHECK =====\n")
cat("Duplicates in old table: ", sum(duplicated(old$Protein_ID)), "\n")
cat("Duplicates in new table: ", sum(duplicated(new$Protein_ID)), "\n\n")

# Check KEGG and COG category columns in new table -----------------------------

cat("===== CATEGORY COLUMN CHECK (NEW TABLE) =====\n")

# crude heuristics to find KEGG and COG category columns
kegg_cols_new <- grep(
  pattern = "^Genetic_|^Metabolism|Transporters|ABC_transporters|Major_facilitator_superfamily__MFS_|Other_transporters|Phosphotransferase_system__PTS_",
  x       = names(new),
  value   = TRUE
)

cog_cols_new <- grep(
  pattern = "^[A-Z]_",
  x       = names(new),
  value   = TRUE
)

cat("Number of KEGG category columns in new: ", length(kegg_cols_new), "\n")
cat("Number of COG category columns in new:  ", length(cog_cols_new), "\n\n")

# Sanity check: logical columns in new should contain only TRUE/FALSE
log_df <- new %>% select(where(is.logical))

logical_ok <- all(sapply(log_df, function(x) all(x %in% c(TRUE, FALSE))))

cat("Logical columns valid (TRUE/FALSE only): ", logical_ok, "\n\n")

cat("===== VALIDATION COMPLETE =====\n")

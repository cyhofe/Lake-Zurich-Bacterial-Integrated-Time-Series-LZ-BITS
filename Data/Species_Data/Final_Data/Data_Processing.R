#clear R's brain
rm(list = ls())

#libraries
library(utils)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)

#set working directory
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Species_Data/")

#load data
Species_Dat <- read.delim(file = "Raw_Data/RepresentativeMAGs.list.txt", header = FALSE, col.names = "MAG_ID")
Checkm2_Dat <- read_tsv(file = "Raw_Data/RepresentativeMagsQualityTable.txt")
GenomicFeatures_Dat <- read.delim(file = "Raw_Data/RepresentativeMagsFeatureTable.txt")
Taxonomy_Dat <- read_tsv(file = "Raw_Data/RepresentativeMAGs_Full_TaxonomyTable.tsv")
rRNA_Dat <- read_delim(file = "Raw_Data/Representative_MAGs.rRNA-SeqReport.txt", col_names = FALSE)
Seqstat_Dat <- read_delim(file = "Raw_Data/RepresentativeMagsSizeTable.Final.txt")
Abundance_Dat <- read_delim(file = "Raw_Data/MergedAbundanceFile.txt")

# prepare the data #################################################################################################

#species dat
Species_Dat <- Species_Dat %>%
  rename("Species_ID" = "MAG_ID")

#checkm2 dat
Checkm2_Dat <- Checkm2_Dat %>%
  rename("Species_ID" = "MAG_ID(1)",
         "Completeness" = "Completeness(2)",
         "Contamination" = "Contamination(3)")

#rRNA dat
rRNA_Dat <- rRNA_Dat %>%
  rename("Species_ID" = "X1",
         "Nr_5s" = "X2",
         "Nr_16s" = "X3",
         "Nr_23s" = "X4")

#Seqstat dat
Seqstat_Dat <- Seqstat_Dat %>%
  rename("Species_ID" = "MAG_ID(1)",
         "Nr_Contigs" = "Nr_Fragments(4)",
         "Nr_Bases" = "MAG_Size_(5)",
         "Min_Contig_Size" = "Min_Contig_Size(6)",
         "Max_Contig_Size" = "Max_Contig_Size(7)",
         "Avg_Contig_Size" = "Avg_Contig_Size(8)")

#Genomic Features

#Select columns
GenomicFeatures_Dat <- GenomicFeatures_Dat %>%
  select(-NUM_CONTIGS, -NUM_5S, -NUM_16S, -NUM_23S, -GENOME_LENGTH)

#remove .gbk from the FILENAME column
GenomicFeatures_Dat <- GenomicFeatures_Dat %>%
  mutate(FILENAME = str_remove(pattern = ".gbk", FILENAME))

# rename columns
GenomicFeatures_Dat <- GenomicFeatures_Dat %>%
  rename(
    "Species_ID" = "FILENAME",
    "GC_Content" = "GC",
    "Total_CDS" = "TOTAL_CDS",
    "Median_CDS" = "MEDIAN_CDS",
    "Mean_CDS" = "MEAN_CDS",
    "Total_tRNA" = "TOTAL_TRNA",
    "Median_tRNA" = "MEDIAN_TRNA",
    "Mean_tRNA" = "MEAN_TRNA",
    "Total_rRNA" = "TOTAL_RRNA",
    "Coding_Density" = "CODING_DENSITY",
    "Noncoding_Density" = "NONCODING_DENSITY",
    "Median_Intergenic_Spacer" = "MEDIAN_INTERGENIC_SPACER",
    "Mean_Intergenic_Spacer" = "MEAN_INTERGENIC_SPACER",
    "Min_Dist" = "MIN_DIST",
    "Max_Dist" = "MAX_DIST",
    "Nr_Features_Used" = "NUM_FEATURES_USED",
    "Negoverlaps_Perc" = "NEGOVERLAPS_PERC",
    "Zerooverlaps_Perc" = "ZEROOVERLAPS_PERC",
    "Posoverlaps_Perc" = "POSOVERLAPS_PERC",
    "Pcc50" = "PCC50",
    "Pcc90" = "PCC90",
    "Nr_Sigma" = "NUM_SIGMA",
    "Nr_Rhodopsins" = "NUM_RHODOPSINS",
    "Nr_Signaltrans" = "NUM_SIGNALTRANS",
    "Hiska" = "HISKA",
    "Ggdef" = "GGDEF",
    "Pas" = "PAS",
    "Pp2c" = "PP2C",
    "Stop_Tag" = "STOP.TAG",
    "Stop_Taa" = "STOP.TAA",
    "Stop_Tga" = "STOP.TGA",
    "Stop_Other" = "STOP.OTHER",
    "Defense_Count" = "DEFENSE_COUNT",
    "Defense_Perc" = "DEFENSE_PERC",
    "Cog_Signal_Trans_Count" = "COG_SIGNAL_TRANS_COUNT",
    "Cog_Signal_Trans_Perc" = "COG_SIGNAL_TRANS_PERC",
    "Phage_Count" = "PHAGE_COUNT",
    "Phage_Perc" = "PHAGE_PERC",
    "Tnp_Count" = "TNP_COUNT",
    "Tnp_Perc" = "TNP_PERC"
  )

GenomicFeatures_Dat <- GenomicFeatures_Dat %>%
  mutate(Total_rRNA = ifelse(is.na(Total_rRNA), 0, Total_rRNA))

# Abundance data

Abundance_Dat <- Abundance_Dat %>%
  pivot_longer(
    cols = -`MAG_ID(1)`,          # Keep MAG_ID(1) as a key column
    names_to = "Date",            # New column for the column names (dates)
    values_to = "Abundance"       # New column for the values (abundances)
  ) %>%
  mutate(
    Date = str_extract(Date, "^[^_]+"),
    Date = as.Date(Date, format = "%d%b%y")) %>%
  rename("Species_ID" = "MAG_ID(1)")

# Taxonomy data

#select important columns
Taxonomy_Dat <- Taxonomy_Dat %>%
  select(user_genome, classification, closest_genome_ani, closest_genome_reference)

#prepare classification
Taxonomy_Dat <- Taxonomy_Dat %>%
  mutate(
    Domain = sub(".*d__([^;]+).*", "\\1", classification),
    Phylum = sub(".*p__([^;]+).*", "\\1", classification),
    Class = sub(".*c__([^;]+).*", "\\1", classification),
    Order = sub(".*o__([^;]+).*", "\\1", classification),
    Family = sub(".*f__([^;]+).*", "\\1", classification),
    Genus = sub(".*g__([^;]+).*", "\\1", classification),
    Species = sub(".*s__([^;]*).*", "\\1", classification)
  ) %>%
  mutate(across(Domain:Species, ~ ifelse(. == classification | . == "", "Unclassified", .))) %>%
  select(-classification)

#put real na's for N/A in column ani
Taxonomy_Dat <- Taxonomy_Dat %>%
  mutate(closest_genome_ani = na_if(closest_genome_ani, "N/A"),
         closest_genome_reference = na_if(closest_genome_reference, "N/A"))

#Rename columns
Taxonomy_Dat <- Taxonomy_Dat %>%
  rename("Species_ID" = "user_genome",
         "ANI_To_Reference" = "closest_genome_ani",
         "Accession_Number" = "closest_genome_reference")

# Count rows with 'Unclassified' in any taxonomic level
unclassified_rows <- Taxonomy_Dat %>%
  filter(if_any(Domain:Species, ~ . == "Unclassified"))

# Add unique tags for unclassified species using last known classification
Taxonomy_Dat <- Taxonomy_Dat %>%
  mutate(
    Genus = ifelse(Genus == "Unclassified" & Family != "Unclassified", paste0(Family, "_unclassified_genus1"), Genus),
    Species = case_when(
      Species == "Unclassified" & Genus != "Unclassified" ~ paste0(Genus, "_unclassified_sp", row_number()),
      Species == "Unclassified" & Family != "Unclassified" ~ paste0(Family, "_unclassified_genus1_sp", row_number()),
      TRUE ~ Species
    )
  )

# Count rows with 'Unclassified' in any taxonomic level
unclassified_rows <- Taxonomy_Dat %>%
  filter(if_any(Domain:Species, ~ . == "Unclassified"))

# Combine the data ##############################################################################################

Dat <- Species_Dat %>%
  left_join(Checkm2_Dat, by = "Species_ID") %>%
  left_join(Seqstat_Dat, by = "Species_ID") %>%
  left_join(Taxonomy_Dat, by = "Species_ID") %>%
  left_join(rRNA_Dat, by = "Species_ID") %>%
  left_join(GenomicFeatures_Dat, by = "Species_ID") %>%
  left_join(Abundance_Dat, by = "Species_ID")

# Calculate Estimated Genome Size
Dat <- Dat %>%
  mutate(Estimated_Genome_Size = (Nr_Bases / Completeness) * (100 - Contamination)) %>%
  relocate(Estimated_Genome_Size, .after = Contamination) %>%
  relocate(Date, .after = Species_ID) %>%
  relocate(Abundance, .after = Date)


# change any N/A in the dataframe to na recognized by R
Dat <- Dat %>%
  mutate(across(where(is.character), ~ ifelse(. == "N/A", NA, .)))

# make a dataframe with each row containing at least one na -> sanity check
Rows_With_NA_In_Dat <- Dat %>%
  filter(if_any(everything(), is.na))

# Calculate total abundance per date
Dat <- Dat %>%
  group_by(Date) %>%
  mutate(Total_Species_Abundance_Per_Date = sum(Abundance))

# Calculate relative species abundance per date
Dat <- Dat %>%
  group_by(Date) %>%
  mutate(Relative_Species_Abundance_Per_Date = (Abundance/Total_Species_Abundance_Per_Date)) %>%
  relocate(Relative_Species_Abundance_Per_Date, .after = Abundance) %>%
  relocate(Total_Species_Abundance_Per_Date, .after = Relative_Species_Abundance_Per_Date)

# rename abundance
Dat <- Dat %>%
  rename("Species_Abundance" = "Abundance")

# Export the data  ###############################################################################
write.csv(x = Dat, file = "Final_Data/Species_Table.csv", row.names = FALSE)

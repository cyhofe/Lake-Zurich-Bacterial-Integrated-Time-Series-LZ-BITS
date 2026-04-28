#clear R's brain
rm(list = ls())

#libraries
library(utils)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)

#set working directory
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/MAG_Data/")

#load data
MAG_IDs <- read.delim(file = "Raw_Data/Filtered_MAGs.list.txt", header = FALSE, col.names = "MAG_ID")
Checkm2_Dat <- read_tsv(file = "Raw_Data/FilteredMAGs_Checkm2_quality_report.tsv")
GenomicFeatures_Dat <- read.delim(file = "Raw_Data/FilteredMAGs_GenomicFeatures.txt")
TaxonomyArch_Dat <- read_tsv(file = "Raw_Data/FilteredMAGs_gtdbtk_ar53_taxonomy.tsv")
TaxonomyBac_Dat <- read_tsv(file = "Raw_Data/FilteredMAGs_gtdbtk_bac120_taxonomy.tsv")
rRNA_Dat <- read_tsv(file = "Raw_Data/FilteredMAGs_rRNA_SeqReport.txt")
Seqstat_Dat <- read_delim(file = "Raw_Data/FilteredMAGs_Seqstat.Report.txt")
Species_Dat <- read_delim(file = "Raw_Data/FilteredMAGs_ClusterFile.txt")

# taxonomy data #################################################################################

#combine the taxonomy data
Taxonomy_Dat <- rbind(TaxonomyArch_Dat, TaxonomyBac_Dat)

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
  mutate(closest_genome_ani = na_if(closest_genome_ani, "N/A"))

#Rename columns
Taxonomy_Dat <- Taxonomy_Dat %>%
  rename("MAG_ID" = "user_genome",
         "ANI_To_Reference" = "closest_genome_ani",
         "Accession_Number" = "closest_genome_reference")

# checkm2 data #################################################################################

#select important columns
Checkm2_Dat <- Checkm2_Dat %>%
  select("Name", "Completeness", "Contamination")

#Rename columns
Checkm2_Dat <- Checkm2_Dat %>%
  rename("MAG_ID" = "Name")

#seqstat data ####################################################################################

#rename columns
Seqstat_Dat <- Seqstat_Dat %>%
  rename("MAG_ID" = "Bin:",
         "Nr_Contigs" = "Number_of_sequences:",
         "Nr_Bases" = "Total_#_residues:",
         "Min_Contig_Size" = "Smallest:",
         "Max_Contig_Size" = "Largest:",
         "Avg_Contig_Size" = "Average_length:")


# GenomicFeatures_Dat ####################################################################################

#Select columns
GenomicFeatures_Dat <- GenomicFeatures_Dat %>%
  select(-NUM_CONTIGS, -NUM_5S, -NUM_16S, -NUM_23S, -GENOME_LENGTH)

#remove .gbk from the FILENAME column
GenomicFeatures_Dat <- GenomicFeatures_Dat %>%
  mutate(FILENAME = str_remove(pattern = ".gbk", FILENAME))

# rename columns
GenomicFeatures_Dat <- GenomicFeatures_Dat %>%
  rename(
    "MAG_ID" = "FILENAME",
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


# Species_Dat ####################################################################################

Species_Dat <- Species_Dat %>%
  rename("MAGs_Per_Species" = "Cluster_Size") %>%
  relocate(MAG_ID, .before = everything()) %>%
  mutate(MAG_ID = str_remove(pattern = ".fna", MAG_ID),
         Species_ID = str_remove(pattern = ".fna", Species_ID))

# Combine the data ###############################################################################

Dat <- Species_Dat %>%
  left_join(Checkm2_Dat, by = "MAG_ID") %>%
  left_join(Seqstat_Dat, by = "MAG_ID") %>%
  left_join(Taxonomy_Dat, by = "MAG_ID") %>%
  left_join(rRNA_Dat, by = "MAG_ID") %>%
  left_join(GenomicFeatures_Dat, by = "MAG_ID")

# Calculate Estimated Genome Size
Dat <- Dat %>%
  mutate(Estimated_Genome_Size = (Dat$Nr_Bases/Dat$Completeness) * (100 - Dat$Contamination)) %>%
  relocate(Estimated_Genome_Size, .after = Contamination)

# change any N/A in the dataframe to na recognized by R
Dat <- Dat %>%
  mutate(across(everything(), ~ ifelse(. == "N/A", NA, .)))

# make a dataframe with each row containing at least one na -> sanity check
Rows_With_NA_In_Dat <- Dat %>%
  filter(if_any(everything(), is.na))

# Export the data  ###############################################################################
write.csv(x = Dat, file = "Final_Data/MAG_Table.csv", row.names = FALSE)

#clear R's brain
rm(list = ls())

#libraries
library(utils)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(readxl)
library(lubridate)
library(data.table)

#set working directory
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Datasets/72Samples_CostumDB/")

#load data
Species_Dat <- read.csv("../../Data/Species_Data/Final_Data/Species_Table.csv")
Peptide_Dat <- read.csv("../../Data/Proteomics_Data/Final_Data/72Samples_p34548_SpringBloom2021_95MagSpeciesDB_15JUL24_20240726_Database_Run07JAN25/Peptide_Table_Combined_Filtered.csv")
Protein_Dat <- read.csv("../../Data/Proteomics_Data/Final_Data/72Samples_p34548_SpringBloom2021_95MagSpeciesDB_15JUL24_20240726_Database_Run07JAN25/Protein_Table_Combined_Filtered.csv")
Lake_Dat <- read.csv("../../Data/Lake_Data/Final_Data/SpringBloom2021_Metadata.csv")
Gene_Dat <- read.csv("../../Data/Gene_Data/Final_Data/Gene_Table.csv")
MAG_Dat <- read.csv("../../Data/MAG_Data/Final_Data/MAG_Table.csv")

# Integrated peptide table
#######################################################################################################################################################
# Prepare Integrated Peptide Table
Integrated_Peptide_Table <- Peptide_Dat %>%
  left_join(Gene_Dat %>% select(Protein_ID, Protein_Sequence) %>% distinct(), by = "Protein_ID") %>%
  left_join(Protein_Dat %>% select(Date, Protein_ID, Protein_Presence, Protein_Group_Size, Nr_Tryptic_Peptides, Protein_Abundance_IBAQ, Relative_Protein_Abundance_IBAQ_Per_Date, Total_Protein_Abundance_IBAQ_Per_Date, Protein_Abundance_MedPolish, Relative_Protein_Abundance_MedPolish_Per_Date, Total_Protein_Abundance_MedPolish_Per_Date) %>% distinct(), by = c("Protein_ID", "Date")) %>%
  left_join(Species_Dat %>% select(Date, Species_ID, Species_Abundance, Relative_Species_Abundance_Per_Date, Total_Species_Abundance_Per_Date, Phylum, Class, Order, Family, Genus, Species), by = c("Date", "Species_ID")) %>%
  left_join(Lake_Dat %>% select(Date, Bacteria_CellsPerMl_Epilimnion) %>% distinct(), by = "Date") %>%
  left_join(MAG_Dat %>% filter(MAG_ID == Species_ID) %>% select(Species_ID, Contamination, Completeness, Estimated_Genome_Size) %>% distinct(), by = "Species_ID") %>%
  distinct()

#manage contaminant entries
Integrated_Peptide_Table <- Integrated_Peptide_Table %>%
  mutate(Phylum = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Phylum),
         Class = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Class),
         Order = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Order),
         Family = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Family),
         Genus = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Genus),
         Species = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), Protein_ID, Species),
         Species_ID = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), Protein_ID, Species_ID))

# sanity check
Summary_Integrated_Peptide_Table <- Integrated_Peptide_Table %>%
  summarize(Nr_Unique_Peptides = n_distinct(Peptide_ID),
            Nr_Unique_Proteins = n_distinct(Protein_ID))

# export the data
write_csv(Integrated_Peptide_Table, "Integrated_Peptide_Table.csv")
#######################################################################################################################################################

# Integrated protein table
#######################################################################################################################################################
# Prepare Integrated Protein Table
Integrated_Protein_Table <- Protein_Dat %>%
  left_join(Gene_Dat %>% select(Protein_ID, KEGG_ID, KEGG_Annotation, COG_ID, COG_Annotation, TIGR_ID, TIGR_Annotation, PFAM_ID, PFAM_Annotation, PHOBIUS_Annotation) %>% distinct(), by = "Protein_ID") %>%
  left_join(Species_Dat %>% select(Date, Species_ID, Species_Abundance, Relative_Species_Abundance_Per_Date, Total_Species_Abundance_Per_Date, Phylum, Class, Order, Family, Genus, Species), by = c("Date", "Species_ID")) %>%
  left_join(Lake_Dat %>% select(Date, Bacteria_CellsPerMl_Epilimnion) %>% distinct(), by = "Date") %>%
  left_join(MAG_Dat %>% filter(MAG_ID == Species_ID) %>% select(Species_ID, Contamination, Completeness, Estimated_Genome_Size) %>% distinct(), by = "Species_ID") %>%
  distinct()

#manage contaminant entries
Integrated_Protein_Table <- Integrated_Protein_Table %>%
  mutate(Phylum = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Phylum),
         Class = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Class),
         Order = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Order),
         Family = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Family),
         Genus = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Genus),
         Species = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), Protein_ID, Species),
         Species_ID = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), Protein_ID, Species_ID))

# sanity check
Summary_Integrated_Protein_Table_1 <- Integrated_Protein_Table %>%
  select(Protein_ID, Nr_Peptides) %>%
  distinct() %>%
  summarize(Nr_Unique_Peptides = sum(Nr_Peptides))

Summary_Integrated_Protein_Table_2 <- Integrated_Protein_Table %>%
  summarize(Nur_Unique_Proteins = n_distinct(Protein_ID))

Summary_Integrated_Protein_Table <- cbind(Summary_Integrated_Protein_Table_1, Summary_Integrated_Protein_Table_2)
rm(Summary_Integrated_Protein_Table_1, Summary_Integrated_Protein_Table_2)

# export the data
write_csv(Integrated_Protein_Table, "Integrated_Protein_Table.csv")
#######################################################################################################################################################

# Integrated genome table
#######################################################################################################################################################
# Prepare Integrated Genome Table
Integrated_Gene_Table <- Gene_Dat %>%
  left_join(Species_Dat %>% select(Species_ID, Phylum, Class, Order, Family, Genus, Species) %>% distinct(), by = c("Species_ID")) %>%
  left_join(MAG_Dat %>% filter(MAG_ID == Species_ID) %>% select(Species_ID, Contamination, Completeness, Estimated_Genome_Size) %>% distinct(), by = "Species_ID") %>%
  distinct()

#manage contaminant entries
Integrated_Gene_Table <- Integrated_Gene_Table %>%
  mutate(Phylum = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Phylum),
         Class = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Class),
         Order = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Order),
         Family = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Family),
         Genus = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), "Y-FGCZCont", Genus),
         Species = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), Protein_ID, Species),
         Species_ID = ifelse(str_detect(Protein_ID, "Y-FGCZCont"), Protein_ID, Species_ID))

# Sanity check
Summary_Integrated_Gene_Table_2 <- Integrated_Gene_Table %>%
  summarize(Nur_Unique_Genes = n_distinct(Protein_ID))

Summary_Integrated_Gene_Table_3 <- Integrated_Gene_Table %>%
  summarize(Nr_Unique_Proteins_WithOutContaminants = n_distinct(Protein_ID))

Summary_Integrated_Gene_Table <- cbind(Summary_Integrated_Gene_Table_2, Summary_Integrated_Gene_Table_3)
rm(Summary_Integrated_Gene_Table_2, Summary_Integrated_Gene_Table_3)

# export the data
write_csv(Integrated_Gene_Table, "Integrated_Gene_Table.csv")

# Generate Tables for each phyla
#######################################################################################################################################################
# Define the dataframes
dataframes <- list(
  Integrated_Protein_Table = Integrated_Protein_Table,
  Integrated_Gene_Table = Integrated_Gene_Table,
  Integrated_Peptide_Table = Integrated_Peptide_Table
)

# Define the output folders
output_folders <- list(
  Integrated_Protein_Table = "Phylum_Protein_Tables",
  Integrated_Gene_Table = "Phylum_Gene_Tables",
  Integrated_Peptide_Table = "Phylum_Peptide_Tables"
)

# Define the list of phyla
Phyla <- Species_Dat %>%
  select(Phylum) %>%
  distinct() %>%
  pull(Phylum)  # Extract phyla as a vector

# Ensure directories exist
for (folder in output_folders) {
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
}

# Export subsets for each phylum into respective folders
for (phylum in Phyla) {
  for (df_name in names(dataframes)) {
    # Subset the dataframe by phylum
    subset_data <- dataframes[[df_name]] %>%
      filter(Phylum == phylum)  # Adjust if "Phylum" column has a different name in your dataframes
    
    # Create the filename
    file_name <- paste0(phylum, "_", df_name, ".csv")
    
    # Save the subset to the appropriate folder
    file_path <- file.path(output_folders[[df_name]], file_name)
    write_csv(subset_data, file_path)
    
    # Print status
    message("Saved: ", file_path)
  }
}



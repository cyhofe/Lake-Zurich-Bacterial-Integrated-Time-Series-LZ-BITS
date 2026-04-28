#clear R's brain
rm(list = ls())

#libraries
library(utils)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)

#set working directory
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Gene_Data/")

#load data
Ids <- read.delim(file = "Raw_Data/MagContigCds.list.txt")
ProteinSequences <- read.delim(file = "Raw_Data/ProteinTable.txt")
ProteinLength <- read.delim(file = "Raw_Data/ProteinLengthTable.txt")
Kegg <- read.delim(file = "Raw_Data/KeggAnnotation.txt", quote = "")
Cog <- read.delim(file = "Raw_Data/CogAnnotation.txt")
Tigr <- read.delim(file = "Raw_Data/TigrAnnotation.txt")
Pfam <- read.delim(file = "Raw_Data/PfamAnnotation.txt")
Phobius <- read.delim(file = "Raw_Data/PhobiusAnnotation.txt")
KEGG_Brite_Categories <- read.csv("Raw_Data/KEGG_Categories/FullCategories/KEGG_Categories_Combined_LongFormat.csv")

# prepare the data #####################################################################################

# prepare transporters and brite categories
Transporters <- KEGG_Brite_Categories %>% 
  ## keep only transporters and the two columns you need
  filter(FunctionalCategory == "Transporters") %>% 
  select(KEGG_ID, FunctionalSubcategory) %>% 
  distinct() %>%                     # each (ID, subcategory) only once
  mutate(value = TRUE) %>%           # flag that will be pivoted
  pivot_wider(
    names_from  = FunctionalSubcategory,
    values_from = value,
    values_fill = FALSE              # what to write where the pair is absent
  )

Brite_Categories <- KEGG_Brite_Categories %>% 
  select(KEGG_ID, FunctionalCategory) %>% 
  distinct() %>% 
  mutate(value = TRUE) %>%                         # flag to widen
  pivot_wider(
    names_from  = FunctionalCategory,
    values_from = value,
    values_fill = FALSE                            # every missing combo → FALSE
  )

#Ids
Ids <- Ids %>%
  rename("Species_ID" = "MAG_ID.1.",
         "Contig_ID" = "Contig_ID.2.",
         "Protein_ID" = "Protein_ID.3.") %>%
  group_by(Species_ID) %>%
  mutate(Nr_Contigs = n_distinct(Contig_ID),
         Nr_Genes = n_distinct(Protein_ID))

#ProteinSequences 
ProteinSequences <- ProteinSequences %>%
  rename("Protein_ID" = "Protein_ID.3.",
         "Protein_Sequence" = "Protein_Sequence.5.") %>%
  mutate(Protein_ID = str_remove(pattern = ">", Protein_ID))

#Proteinlength
ProteinLength <- ProteinLength  %>%
  rename("Protein_ID" = "Protein_ID.3.",
         "Protein_Length" = "Protein_Length.6.") %>%
  mutate(Protein_ID = str_remove(pattern = ">", Protein_ID))

#Kegg

# rename
Kegg <- Kegg %>%
  rename("Protein_ID" = "Protein_ID.3.",
         "KEGG_Annotation" = "KEGG_Annotation.7.")

# correct annotations
Kegg <- Kegg %>%
  mutate(
    KEGG_Annotation = str_replace(pattern = "-:-", replacement = "Unclassified", KEGG_Annotation),
    KEGG_ID = str_extract(KEGG_Annotation, "^[^:]+"),  # Extract everything before the first ":"
    KEGG_Annotation = str_extract(KEGG_Annotation, "(?<=:).*"),
    KEGG_Annotation = ifelse(is.na(KEGG_Annotation), "UnknownFunction", KEGG_Annotation))



# add kegg categories
Kegg <- Kegg %>%
  left_join(Brite_Categories, by = "KEGG_ID") %>%
  left_join(Transporters, by = "KEGG_ID")

#Cog

#rename 
Cog <- Cog %>%
  rename("Protein_ID" = "Protein_ID.3.",
         "COG_Annotation" = "COG_Annotation.9.")

# correct annotations
Cog <- Cog %>%
  mutate(
    COG_Annotation = str_replace(pattern = "-_-", replacement = "Unclassified", COG_Annotation),
    COG_ID = str_extract(COG_Annotation, "^[^:]+"),  # Extract everything before the first ":"
    COG_Annotation = str_extract(COG_Annotation, "(?<=:).*"),
    COG_Annotation = ifelse(is.na(COG_Annotation), "UnknownFunction", COG_Annotation))

#Tigr

#rename 
Tigr <- Tigr %>%
  rename("Protein_ID" = "Protein_ID.3.",
         "TIGR_Annotation" = "TIGR_Annotation.8.")

# correct annotations
Tigr <- Tigr %>%
  mutate(
    TIGR_Annotation = str_replace(pattern = "-_-", replacement = "Unclassified", TIGR_Annotation),
    TIGR_ID = str_extract(TIGR_Annotation, "^[^:]+"),  # Extract everything before the first ":"
    TIGR_Annotation = str_extract(TIGR_Annotation, "(?<=:).*"),
    TIGR_Annotation = ifelse(is.na(TIGR_Annotation), "UnknownFunction", TIGR_Annotation))

#Pfam

#rename 
Pfam <- Pfam %>%
  rename("Protein_ID" = "Protein_ID.3.",
         "PFAM_Annotation" = "PFAM_Annotation.10.")

# correct annotations
Pfam <- Pfam %>%
  mutate(
    PFAM_Annotation = str_replace(pattern = "-:-", replacement = "Unclassified", PFAM_Annotation),
    PFAM_ID = str_extract(PFAM_Annotation, "^[^:]+"),  # Extract everything before the first ":"
    PFAM_Annotation = str_extract(PFAM_Annotation, "(?<=:).*"),
    PFAM_Annotation = ifelse(is.na(PFAM_Annotation), "UnknownFunction", PFAM_Annotation))


# Phobius
Phobius <- Phobius %>%
  rename("Protein_ID" = "Protein_ID.3.",
         "PHOBIUS_Annotation" = "PHOBIUS_Annotation.11.")

# correct annotations
Phobius <- Phobius %>%
  mutate(PHOBIUS_Annotation = ifelse(PHOBIUS_Annotation == "Y", "SecretedProtein", PHOBIUS_Annotation),
         PHOBIUS_Annotation = ifelse(PHOBIUS_Annotation == "0", "CytoplasmaticProtein", PHOBIUS_Annotation))


################### combine the data ####################################################################################

Dat <- Ids %>%
  left_join(Kegg, by = "Protein_ID") %>%
  left_join(Cog, by = "Protein_ID") %>%
  left_join(Tigr, by = "Protein_ID") %>%
  left_join(Pfam, by = "Protein_ID") %>%
  left_join(Phobius, by = "Protein_ID") %>%
  left_join(ProteinLength, by = "Protein_ID") %>%
  left_join(ProteinSequences, by = "Protein_ID") %>%
  select(-c(8:35), c(8:35))

################### Fuse Contig and Protein IDs ########################################################################

# fix protein names according to database
Dat <- Dat %>%
  mutate(
    Protein_ID = paste0(
      Species_ID,                               # Add Species_ID
      "-",                                      # Add separator
      str_extract(Contig_ID, "contig_\\d+"),    # Extract contig number
      "-",                                      # Add separator
      str_extract(Protein_ID, "cds\\d+")        # Extract cds number
    )
  )

# replace missing values with False to mention that these dont belong to these categories
Dat <- Dat %>%
  mutate(across(
    c(
      Genetic_Information_Processing,
      Metabolism,
      Signaling_and_Cellular_Processes,
      Transporters,
      `ABC_transporters__prokaryotic_type - ABC-2_type_and_other_transporters`,
      `ABC_transporters__prokaryotic_type - Metallic_cation__iron-siderophore_and_vitamin_B12_transporters`,
      `ABC_transporters__prokaryotic_type - Mineral_and_organic_ion_transporters`,
      `ABC_transporters__prokaryotic_type - Peptide_and_nickel_transporters`,
      `ABC_transporters__prokaryotic_type - Phosphate_and_amino_acid_transporters`,
      `ABC_transporters__prokaryotic_type - Saccharide__polyol__and_lipid_transporters`,
      `Major_facilitator_superfamily__MFS_ - Drug_transporters`,
      `Major_facilitator_superfamily__MFS_ - Lipid_transporters`,
      `Major_facilitator_superfamily__MFS_ - Metal_transporters`,
      `Major_facilitator_superfamily__MFS_ - Nitrate_nitrite_transporters`,
      `Major_facilitator_superfamily__MFS_ - Organic_acid_transporters`,
      `Major_facilitator_superfamily__MFS_ - Phosphate_and_organophosphate_transporters`,
      `Major_facilitator_superfamily__MFS_ - Protein_transporters`,
      `Major_facilitator_superfamily__MFS_ - Sugar_transporters`,
      `Major_facilitator_superfamily__MFS_ - Unknown_transporters`,
      `Other_transporters - Accessory_factors_involved_in_transport_[TC_8]`,
      `Other_transporters - Aquaporins_and_small_neutral_solute_transporters_[TC_1_A_8]`,
      `Other_transporters - Electrochemical_potential-driven_transporters_[TC_2]`,
      `Other_transporters - Others`,
      `Other_transporters - Pores_ion_channels_[TC_1]`,
      `Other_transporters - Primary_active_transporters_[TC_3]`,
      `Other_transporters - Transmembrane_electron_carriers_[TC_5]`,
      `Phosphotransferase_system__PTS_ - Enzyme_I_and_HPr`,
      `Phosphotransferase_system__PTS_ - Enzyme_II_[TC_4_A]`
    ),
    ~ replace_na(.x, FALSE)
  ))



# export the data
write.csv(Dat,file = "Final_Data/Gene_Table.csv", row.names = FALSE)


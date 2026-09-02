# clear R's brain
rm(list = ls())

# load libraries
library(dplyr)
library(stringr)
library(readr)

# set working directory
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Method_Paper_Analysis/FinalAnalysis/Figure_03/MockCommunityAnalysis/Database/")

# read data
raw_dat <- read.csv(file = "MockD6300.Genes.csv", header = F)

# prepare data
dat <- raw_dat %>%
  select(V1, V2) %>%
  rename("Gene_ID" = "V1",
         "Annotation" = "V2")

# prepare species id
species <- dat %>%
  select(Gene_ID) %>%
  mutate(Species = case_when(str_detect(Gene_ID, "Cont") ~ "Contaminant",
                             str_detect(Gene_ID, "fgcz") ~ "Contaminant",
                             TRUE ~ str_extract(Gene_ID, "^[^_]+_[^_]+")))

# export data
readr::write_csv(species, "MockD6300.Genes.Perpared.csv")

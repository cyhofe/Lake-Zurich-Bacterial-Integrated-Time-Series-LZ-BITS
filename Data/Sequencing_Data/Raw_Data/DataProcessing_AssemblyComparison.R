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
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Sequencing_Data/Raw_Data/")

# Load Illumina Raw ##################################################################################
#load data
folder_path <- "Assembly_Comparison_Statistics/Illumina_Raw/"  

# ——— find all .txt files ———
file_paths <- list.files(
  path       = folder_path,
  pattern    = "\\.txt$",
  full.names = TRUE
)

# ——— a helper to parse one stats‐file ———
parse_stats_file <- function(fp) {
  # read as two columns: key and value+comments
  kv <- read.table(fp,
                   sep            = "\t",
                   header         = FALSE,
                   stringsAsFactors = FALSE,
                   fill           = TRUE,
                   comment.char   = "")  
  
  # clean up key and numeric value
  keys   <- gsub(":$", "", kv[[1]])                          # strip trailing colon
  vals   <- as.numeric(gsub("\\s*#.*$", "", kv[[2]]))        # drop inline “# …” then coerce
  
  stats  <- setNames(vals, keys)
  
  # pull out what you asked for
  bases_mapped   <- stats["bases mapped"]
  total_reads    <- stats["sequences"]
  reads_mapped   <- stats["reads mapped"]
  mismatches     <- stats["mismatches"]
  error_rate     <- stats["error rate"]
  
  # compute percent mapped
  pct_mapped     <- reads_mapped / total_reads * 100
  
  # return a one‐row data.frame
  data.frame(
    sample         = tools::file_path_sans_ext(basename(fp)),
    total_reads    = total_reads,
    reads_mapped   = reads_mapped,
    percent_mapped = pct_mapped,
    mismatches     = mismatches,
    error_rate     = error_rate,
    bases_mapped   = bases_mapped,
    stringsAsFactors = FALSE)}

# ——— apply to all files and bind into one table ———
stats_list <- lapply(file_paths, parse_stats_file)
summary_df <- do.call(rbind, stats_list)

# modify and add sample and date info
Illumina_Raw <- summary_df %>%
  mutate(Sample = word(sample, 1, sep = fixed(".")),
         Assembly = "Illumina_Raw",
         Date = as.Date(str_extract(Sample, "[0-9]{2}[A-Za-z]{3}[0-9]{2}"), format = "%d%b%y")) %>% 
  select(-sample)


# Load Illumina Filtered ##################################################################################
#load data
folder_path <- "Assembly_Comparison_Statistics/Illumina_Filtered/"  

# ——— find all .txt files ———
file_paths <- list.files(
  path       = folder_path,
  pattern    = "\\.txt$",
  full.names = TRUE
)

# ——— a helper to parse one stats‐file ———
parse_stats_file <- function(fp) {
  # read as two columns: key and value+comments
  kv <- read.table(fp,
                   sep            = "\t",
                   header         = FALSE,
                   stringsAsFactors = FALSE,
                   fill           = TRUE,
                   comment.char   = "")  
  
  # clean up key and numeric value
  keys   <- gsub(":$", "", kv[[1]])                          # strip trailing colon
  vals   <- as.numeric(gsub("\\s*#.*$", "", kv[[2]]))        # drop inline “# …” then coerce
  
  stats  <- setNames(vals, keys)
  
  # pull out what you asked for
  bases_mapped   <- stats["bases mapped"]
  total_reads    <- stats["sequences"]
  reads_mapped   <- stats["reads mapped"]
  mismatches     <- stats["mismatches"]
  error_rate     <- stats["error rate"]
  
  # compute percent mapped
  pct_mapped     <- reads_mapped / total_reads * 100
  
  # return a one‐row data.frame
  data.frame(
    sample         = tools::file_path_sans_ext(basename(fp)),
    total_reads    = total_reads,
    reads_mapped   = reads_mapped,
    percent_mapped = pct_mapped,
    mismatches     = mismatches,
    error_rate     = error_rate,
    bases_mapped   = bases_mapped,
    stringsAsFactors = FALSE)}

# ——— apply to all files and bind into one table ———
stats_list <- lapply(file_paths, parse_stats_file)
summary_df <- do.call(rbind, stats_list)

# modify and add sample and date info
Illumina_Filtered <- summary_df %>%
  mutate(Sample = word(sample, 1, sep = fixed(".")),
         Assembly = "Illumina_Filtered",
         Date = as.Date(str_extract(Sample, "[0-9]{2}[A-Za-z]{3}[0-9]{2}"), format = "%d%b%y")) %>% 
  select(-sample)

# Load Nanopore Raw ##################################################################################
#load data
folder_path <- "Assembly_Comparison_Statistics/Nanopore_Raw/"  

# ——— find all .txt files ———
file_paths <- list.files(
  path       = folder_path,
  pattern    = "\\.txt$",
  full.names = TRUE
)

# ——— a helper to parse one stats‐file ———
parse_stats_file <- function(fp) {
  # read as two columns: key and value+comments
  kv <- read.table(fp,
                   sep            = "\t",
                   header         = FALSE,
                   stringsAsFactors = FALSE,
                   fill           = TRUE,
                   comment.char   = "")  
  
  # clean up key and numeric value
  keys   <- gsub(":$", "", kv[[1]])                          # strip trailing colon
  vals   <- as.numeric(gsub("\\s*#.*$", "", kv[[2]]))        # drop inline “# …” then coerce
  
  stats  <- setNames(vals, keys)
  
  # pull out what you asked for
  bases_mapped   <- stats["bases mapped"]
  total_reads    <- stats["sequences"]
  reads_mapped   <- stats["reads mapped"]
  mismatches     <- stats["mismatches"]
  error_rate     <- stats["error rate"]
  
  # compute percent mapped
  pct_mapped     <- reads_mapped / total_reads * 100
  
  # return a one‐row data.frame
  data.frame(
    sample         = tools::file_path_sans_ext(basename(fp)),
    total_reads    = total_reads,
    reads_mapped   = reads_mapped,
    percent_mapped = pct_mapped,
    mismatches     = mismatches,
    error_rate     = error_rate,
    bases_mapped   = bases_mapped,
    stringsAsFactors = FALSE)}

# ——— apply to all files and bind into one table ———
stats_list <- lapply(file_paths, parse_stats_file)
summary_df <- do.call(rbind, stats_list)

# modify and add sample and date info
Nanopore_Raw <- summary_df %>%
  mutate(Sample = word(sample, 1, sep = fixed(".")),
         Assembly = "Nanopore_Raw",
         Date = as.Date(str_extract(Sample, "[0-9]{2}[A-Za-z]{3}[0-9]{2}"), format = "%d%b%y")) %>% 
  select(-sample)

# Load Nanopore Filtered ##################################################################################
#load data
folder_path <- "Assembly_Comparison_Statistics/Nanopore_Filtered/"  

# ——— find all .txt files ———
file_paths <- list.files(
  path       = folder_path,
  pattern    = "\\.txt$",
  full.names = TRUE
)

# ——— a helper to parse one stats‐file ———
parse_stats_file <- function(fp) {
  # read as two columns: key and value+comments
  kv <- read.table(fp,
                   sep            = "\t",
                   header         = FALSE,
                   stringsAsFactors = FALSE,
                   fill           = TRUE,
                   comment.char   = "")  
  
  # clean up key and numeric value
  keys   <- gsub(":$", "", kv[[1]])                          # strip trailing colon
  vals   <- as.numeric(gsub("\\s*#.*$", "", kv[[2]]))        # drop inline “# …” then coerce
  
  stats  <- setNames(vals, keys)
  
  # pull out what you asked for
  bases_mapped   <- stats["bases mapped"]
  total_reads    <- stats["sequences"]
  reads_mapped   <- stats["reads mapped"]
  mismatches     <- stats["mismatches"]
  error_rate     <- stats["error rate"]
  
  # compute percent mapped
  pct_mapped     <- reads_mapped / total_reads * 100
  
  # return a one‐row data.frame
  data.frame(
    sample         = tools::file_path_sans_ext(basename(fp)),
    total_reads    = total_reads,
    reads_mapped   = reads_mapped,
    percent_mapped = pct_mapped,
    mismatches     = mismatches,
    error_rate     = error_rate,
    bases_mapped   = bases_mapped,
    stringsAsFactors = FALSE)}

# ——— apply to all files and bind into one table ———
stats_list <- lapply(file_paths, parse_stats_file)
summary_df <- do.call(rbind, stats_list)

# modify and add sample and date info
Nanopore_Filtered <- summary_df %>%
  mutate(Sample = word(sample, 1, sep = fixed(".")),
         Assembly = "Nanopore_Filtered",
         Date = as.Date(str_extract(Sample, "[0-9]{2}[A-Za-z]{3}[0-9]{2}"), format = "%d%b%y")) %>% 
  select(-sample)


# Load Nanopore Corrected ##################################################################################
#load data
folder_path <- "Assembly_Comparison_Statistics/Nanopore_Corrected/"  

# ——— find all .txt files ———
file_paths <- list.files(
  path       = folder_path,
  pattern    = "\\.txt$",
  full.names = TRUE
)

# ——— a helper to parse one stats‐file ———
parse_stats_file <- function(fp) {
  # read as two columns: key and value+comments
  kv <- read.table(fp,
                   sep            = "\t",
                   header         = FALSE,
                   stringsAsFactors = FALSE,
                   fill           = TRUE,
                   comment.char   = "")  
  
  # clean up key and numeric value
  keys   <- gsub(":$", "", kv[[1]])                          # strip trailing colon
  vals   <- as.numeric(gsub("\\s*#.*$", "", kv[[2]]))        # drop inline “# …” then coerce
  
  stats  <- setNames(vals, keys)
  
  # pull out what you asked for
  bases_mapped   <- stats["bases mapped"]
  total_reads    <- stats["sequences"]
  reads_mapped   <- stats["reads mapped"]
  mismatches     <- stats["mismatches"]
  error_rate     <- stats["error rate"]
  
  # compute percent mapped
  pct_mapped     <- reads_mapped / total_reads * 100
  
  # return a one‐row data.frame
  data.frame(
    sample         = tools::file_path_sans_ext(basename(fp)),
    total_reads    = total_reads,
    reads_mapped   = reads_mapped,
    percent_mapped = pct_mapped,
    mismatches     = mismatches,
    error_rate     = error_rate,
    bases_mapped   = bases_mapped,
    stringsAsFactors = FALSE)}

# ——— apply to all files and bind into one table ———
stats_list <- lapply(file_paths, parse_stats_file)
summary_df <- do.call(rbind, stats_list)

# modify and add sample and date info
Nanopore_Corrected <- summary_df %>%
  mutate(Sample = word(sample, 1, sep = fixed(".")),
         Assembly = "Nanopore_Corrected",
         Date = as.Date(str_extract(Sample, "[0-9]{2}[A-Za-z]{3}[0-9]{2}"), format = "%d%b%y")) %>% 
  select(-sample)



##########################################################
##########################################################

# ——— Illumina_Raw NanoStat files ———
folder_path <- "Assembly_Comparison_Statistics/Illumina_Raw_ContigsInfo/"

# find all .txt (NanoStat) files
file_paths <- list.files(
  path       = folder_path,
  pattern    = "\\.txt$",
  full.names = TRUE
)

# helper to parse one NanoStat report
parse_nanostat <- function(fp) {
  # read two-column table with header “Metrics” / “dataset”
  df <- read.table(fp,
                   sep             = "\t",
                   header          = TRUE,
                   stringsAsFactors = FALSE)
  metrics <- setNames(df$dataset, df$Metrics)
  
  # pull out the 5 fields you care about
  n_reads   <- as.numeric(metrics["number_of_reads"])
  n_bases   <- as.numeric(metrics["number_of_bases"])
  med_len   <- as.numeric(metrics["median_read_length"])
  mean_len  <- as.numeric(metrics["mean_read_length"])
  stdev_len <- as.numeric(metrics["read_length_stdev"])
  
  # extract Sample and Date from the filename
  base   <- tools::file_path_sans_ext(basename(fp))
  sample <- word(base, 1, sep = fixed("."))
  date   <- as.Date(
    str_extract(sample, "[0-9]{2}[A-Za-z]{3}[0-9]{2}"),
    format = "%d%b%y"
  )
  
  tibble(
    Sample             = sample,
    Date               = date,
    number_of_reads    = n_reads,
    number_of_bases    = n_bases,
    median_read_length = med_len,
    mean_read_length   = mean_len,
    read_length_stdev  = stdev_len
  )
}

# apply to all files and combine
Illumina_Raw_NanoStat <- file_paths %>%
  lapply(parse_nanostat) %>%
  bind_rows() %>%
  mutate(Assembly = "Illumina_Raw") %>%
  select(Sample, Date, Assembly,
         number_of_reads, number_of_bases,
         median_read_length, mean_read_length, read_length_stdev)

# peek
print(Illumina_Raw_NanoStat)


# ——— Illumina_Filtered NanoStat files ———
folder_path <- "Assembly_Comparison_Statistics/Illumina_Filtered_ContigsInfo/"

# find all .txt (NanoStat) files
file_paths <- list.files(
  path       = folder_path,
  pattern    = "\\.txt$",
  full.names = TRUE
)

# helper to parse one NanoStat report
parse_nanostat <- function(fp) {
  # read two-column table with header “Metrics” / “dataset”
  df <- read.table(fp,
                   sep             = "\t",
                   header          = TRUE,
                   stringsAsFactors = FALSE)
  metrics <- setNames(df$dataset, df$Metrics)
  
  # pull out the 5 fields you care about
  n_reads   <- as.numeric(metrics["number_of_reads"])
  n_bases   <- as.numeric(metrics["number_of_bases"])
  med_len   <- as.numeric(metrics["median_read_length"])
  mean_len  <- as.numeric(metrics["mean_read_length"])
  stdev_len <- as.numeric(metrics["read_length_stdev"])
  
  # extract Sample and Date from the filename
  base   <- tools::file_path_sans_ext(basename(fp))
  sample <- word(base, 1, sep = fixed("."))
  date   <- as.Date(
    str_extract(sample, "[0-9]{2}[A-Za-z]{3}[0-9]{2}"),
    format = "%d%b%y"
  )
  
  tibble(
    Sample             = sample,
    Date               = date,
    number_of_reads    = n_reads,
    number_of_bases    = n_bases,
    median_read_length = med_len,
    mean_read_length   = mean_len,
    read_length_stdev  = stdev_len
  )
}

# apply to all files and combine
Illumina_Filtered_NanoStat <- file_paths %>%
  lapply(parse_nanostat) %>%
  bind_rows() %>%
  mutate(Assembly = "Illumina_Filtered") %>%
  select(Sample, Date, Assembly,
         number_of_reads, number_of_bases,
         median_read_length, mean_read_length, read_length_stdev)

# peek
print(Illumina_Filtered_NanoStat)

# ——— Nanopore_Raw NanoStat files ———
folder_path <- "Assembly_Comparison_Statistics/Nanopore_Raw_ContigsInfo/"

# find all .txt (NanoStat) files
file_paths <- list.files(
  path       = folder_path,
  pattern    = "\\.txt$",
  full.names = TRUE
)

# helper to parse one NanoStat report
parse_nanostat <- function(fp) {
  # read two-column table with header “Metrics” / “dataset”
  df <- read.table(fp,
                   sep             = "\t",
                   header          = TRUE,
                   stringsAsFactors = FALSE)
  metrics <- setNames(df$dataset, df$Metrics)
  
  # pull out the 5 fields you care about
  n_reads   <- as.numeric(metrics["number_of_reads"])
  n_bases   <- as.numeric(metrics["number_of_bases"])
  med_len   <- as.numeric(metrics["median_read_length"])
  mean_len  <- as.numeric(metrics["mean_read_length"])
  stdev_len <- as.numeric(metrics["read_length_stdev"])
  
  # extract Sample and Date from the filename
  base   <- tools::file_path_sans_ext(basename(fp))
  sample <- word(base, 1, sep = fixed("."))
  date   <- as.Date(
    str_extract(sample, "[0-9]{2}[A-Za-z]{3}[0-9]{2}"),
    format = "%d%b%y"
  )
  
  tibble(
    Sample             = sample,
    Date               = date,
    number_of_reads    = n_reads,
    number_of_bases    = n_bases,
    median_read_length = med_len,
    mean_read_length   = mean_len,
    read_length_stdev  = stdev_len
  )
}

# apply to all files and combine
Nanopore_Raw_NanoStat <- file_paths %>%
  lapply(parse_nanostat) %>%
  bind_rows() %>%
  mutate(Assembly = "Nanopore_Raw") %>%
  select(Sample, Date, Assembly,
         number_of_reads, number_of_bases,
         median_read_length, mean_read_length, read_length_stdev)

# peek
print(Nanopore_Raw_NanoStat)


# ——— Nanopore_Filtered NanoStat files ———
folder_path <- "Assembly_Comparison_Statistics/Nanopore_Filtered_ContigsInfo/"

# find all .txt (NanoStat) files
file_paths <- list.files(
  path       = folder_path,
  pattern    = "\\.txt$",
  full.names = TRUE
)

# helper to parse one NanoStat report
parse_nanostat <- function(fp) {
  # read two-column table with header “Metrics” / “dataset”
  df <- read.table(fp,
                   sep             = "\t",
                   header          = TRUE,
                   stringsAsFactors = FALSE)
  metrics <- setNames(df$dataset, df$Metrics)
  
  # pull out the 5 fields you care about
  n_reads   <- as.numeric(metrics["number_of_reads"])
  n_bases   <- as.numeric(metrics["number_of_bases"])
  med_len   <- as.numeric(metrics["median_read_length"])
  mean_len  <- as.numeric(metrics["mean_read_length"])
  stdev_len <- as.numeric(metrics["read_length_stdev"])
  
  # extract Sample and Date from the filename
  base   <- tools::file_path_sans_ext(basename(fp))
  sample <- word(base, 1, sep = fixed("."))
  date   <- as.Date(
    str_extract(sample, "[0-9]{2}[A-Za-z]{3}[0-9]{2}"),
    format = "%d%b%y"
  )
  
  tibble(
    Sample             = sample,
    Date               = date,
    number_of_reads    = n_reads,
    number_of_bases    = n_bases,
    median_read_length = med_len,
    mean_read_length   = mean_len,
    read_length_stdev  = stdev_len
  )
}

# apply to all files and combine
Nanopore_Filtered_NanoStat <- file_paths %>%
  lapply(parse_nanostat) %>%
  bind_rows() %>%
  mutate(Assembly = "Nanopore_Filtered") %>%
  select(Sample, Date, Assembly,
         number_of_reads, number_of_bases,
         median_read_length, mean_read_length, read_length_stdev)

# peek
print(Nanopore_Filtered_NanoStat)


# ——— Nanopore_Corrected NanoStat files ———
folder_path <- "Assembly_Comparison_Statistics/Nanopore_Corrected_ContigsInfo/"

# find all .txt (NanoStat) files
file_paths <- list.files(
  path       = folder_path,
  pattern    = "\\.txt$",
  full.names = TRUE
)

# helper to parse one NanoStat report
parse_nanostat <- function(fp) {
  # read two-column table with header “Metrics” / “dataset”
  df <- read.table(fp,
                   sep             = "\t",
                   header          = TRUE,
                   stringsAsFactors = FALSE)
  metrics <- setNames(df$dataset, df$Metrics)
  
  # pull out the 5 fields you care about
  n_reads   <- as.numeric(metrics["number_of_reads"])
  n_bases   <- as.numeric(metrics["number_of_bases"])
  med_len   <- as.numeric(metrics["median_read_length"])
  mean_len  <- as.numeric(metrics["mean_read_length"])
  stdev_len <- as.numeric(metrics["read_length_stdev"])
  
  # extract Sample and Date from the filename
  base   <- tools::file_path_sans_ext(basename(fp))
  sample <- word(base, 1, sep = fixed("."))
  date   <- as.Date(
    str_extract(sample, "[0-9]{2}[A-Za-z]{3}[0-9]{2}"),
    format = "%d%b%y"
  )
  
  tibble(
    Sample             = sample,
    Date               = date,
    number_of_reads    = n_reads,
    number_of_bases    = n_bases,
    median_read_length = med_len,
    mean_read_length   = mean_len,
    read_length_stdev  = stdev_len
  )
}

# apply to all files and combine
Nanopore_Corrected_NanoStat <- file_paths %>%
  lapply(parse_nanostat) %>%
  bind_rows() %>%
  mutate(Assembly = "Nanopore_Corrected") %>%
  select(Sample, Date, Assembly,
         number_of_reads, number_of_bases,
         median_read_length, mean_read_length, read_length_stdev)

# peek
print(Nanopore_Filtered_NanoStat)

########################################## combine to one dataframe
dat <- rbind(Illumina_Raw, Illumina_Filtered, Nanopore_Raw, Nanopore_Filtered, Nanopore_Corrected)
dat_nanostat <- rbind(Illumina_Raw_NanoStat, Illumina_Filtered_NanoStat, Nanopore_Raw_NanoStat, Nanopore_Filtered_NanoStat, Nanopore_Corrected_NanoStat)

final_dat <- dat %>%
  left_join(dat_nanostat, by = c("Sample", "Date", "Assembly"))

#set working directory
setwd("/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Sequencing_Data/Raw_Data/")

# export the data
write_csv(final_dat, "AssemblyComparison.csv")

############################################################
# Quality-filtering statistics:
# FastQC (Illumina) + NanoStat (Nanopore)
############################################################

# Clear R's brain
rm(list = ls())


############################################################
# 1) Libraries
############################################################

library(utils)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)


############################################################
# 2) Set working directory
############################################################

setwd(
  "/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Sequencing_Data/"
)

# If R is running directly inside WSL, use this path instead:
# setwd(
#   "/mnt/c/SynologyDrive/PhD/Projects/SpringBloom2021/Data/Sequencing_Data/"
# )


############################################################
# 3) Input and output directories
############################################################

quality_dir <- file.path(
  "Raw_Data",
  "Quality_Trimming"
)

illumina_raw_dir <- file.path(
  quality_dir,
  "Illumina_Raw_FastQC"
)

illumina_filtered_dir <- file.path(
  quality_dir,
  "Illumina_Filtered_FastQC"
)

nanopore_dir <- file.path(
  quality_dir,
  "Nanopore_QualityReports_NanoStat"
)

output_file <- file.path(
  quality_dir,
  "QualityFilteringStats.csv"
)

# Stop immediately if one of the expected folders is missing
required_directories <- c(
  illumina_raw_dir,
  illumina_filtered_dir,
  nanopore_dir
)

missing_directories <- required_directories[
  !dir.exists(required_directories)
]

if (length(missing_directories) > 0L) {
  stop(
    "The following input directories were not found:\n",
    paste(missing_directories, collapse = "\n")
  )
}


############################################################
# 4) General helper functions
############################################################

# Convert values such as:
# "133,324,894" -> 133324894
# "4,279.5"     -> 4279.5
parse_plain_number <- function(x) {
  suppressWarnings(
    parse_number(
      str_replace_all(as.character(x), ",", ""),
      na = c("", "NA", "NaN")
    )
  )
}


# Convert FastQC total-base values such as:
# "19.9 Gbp" -> 19900000000
# "850 Mbp"  -> 850000000
parse_scaled_bases <- function(x) {
  
  clean_value <- str_to_lower(str_trim(as.character(x)))
  numeric_value <- parse_plain_number(clean_value)
  
  unit <- str_match(
    clean_value,
    "(tbp|gbp|mbp|kbp|bp|tb|gb|mb|kb)\\s*$"
  )[, 2]
  
  multiplier <- case_when(
    unit %in% c("kbp", "kb") ~ 1e3,
    unit %in% c("mbp", "mb") ~ 1e6,
    unit %in% c("gbp", "gb") ~ 1e9,
    unit %in% c("tbp", "tb") ~ 1e12,
    TRUE                       ~ 1
  )
  
  numeric_value * multiplier
}


# Convert FastQC bins to a numeric value.
#
# Examples:
# "150"    -> 150
# "50-54"  -> 52
# "10-14"  -> 12
#
# For interval bins, the midpoint is used.
parse_bin_midpoint <- function(x) {
  
  extracted_numbers <- str_extract_all(
    as.character(x),
    "\\d+(?:\\.\\d+)?"
  )
  
  vapply(
    extracted_numbers,
    FUN = function(values) {
      
      values <- as.numeric(values)
      
      if (length(values) == 0L) {
        return(NA_real_)
      }
      
      if (length(values) == 1L) {
        return(values[[1]])
      }
      
      mean(values[1:2])
    },
    FUN.VALUE = numeric(1)
  )
}


# Weighted mean with explicit handling of missing values
weighted_mean_safe <- function(x, weights) {
  
  keep <- is.finite(x) &
    is.finite(weights) &
    weights > 0
  
  if (!any(keep)) {
    return(NA_real_)
  }
  
  weighted.mean(
    x = x[keep],
    w = weights[keep]
  )
}


# Weighted median for binned FastQC distributions
weighted_median <- function(x, weights) {
  
  keep <- is.finite(x) &
    is.finite(weights) &
    weights > 0
  
  if (!any(keep)) {
    return(NA_real_)
  }
  
  x <- x[keep]
  weights <- weights[keep]
  
  ordering <- order(x)
  
  x <- x[ordering]
  weights <- weights[ordering]
  
  cumulative_weights <- cumsum(weights)
  median_position <- sum(weights) / 2
  
  x[
    which(cumulative_weights >= median_position)[1]
  ]
}


# Standardize column names so that names such as
# "Total Sequences", "total_sequences", and "Total.Sequences"
# can be compared safely.
normalize_column_name <- function(x) {
  str_to_lower(
    str_replace_all(
      x,
      "[^[:alnum:]]",
      ""
    )
  )
}


# Find a column in a parsed table without relying on punctuation
find_table_column <- function(data, requested_name) {
  
  normalized_names <- normalize_column_name(names(data))
  normalized_request <- normalize_column_name(requested_name)
  
  matching_column <- which(
    normalized_names == normalized_request
  )
  
  if (length(matching_column) == 0L) {
    stop(
      "Could not find column '",
      requested_name,
      "'. Available columns are: ",
      paste(names(data), collapse = ", ")
    )
  }
  
  names(data)[matching_column[[1]]]
}


############################################################
# 5) FastQC module parser
############################################################

# Read one complete FastQC module by its module name.
#
# This does not depend on fixed line numbers. Therefore, it
# remains valid when modules contain different numbers of rows.
read_fastqc_module <- function(lines, module_name) {
  
  module_start <- which(
    startsWith(
      lines,
      paste0(">>", module_name)
    )
  )
  
  if (length(module_start) == 0L) {
    stop(
      "FastQC module '",
      module_name,
      "' was not found."
    )
  }
  
  module_start <- module_start[[1]]
  
  module_end <- which(
    seq_along(lines) > module_start &
      startsWith(lines, ">>END_MODULE")
  )
  
  if (length(module_end) == 0L) {
    stop(
      "No END_MODULE marker was found after FastQC module '",
      module_name,
      "'."
    )
  }
  
  module_end <- module_end[[1]]
  
  module_lines <- lines[
    seq.int(
      from = module_start + 1L,
      to = module_end - 1L
    )
  ]
  
  header_index <- which(
    startsWith(module_lines, "#")
  )
  
  if (length(header_index) == 0L) {
    stop(
      "No table header was found in FastQC module '",
      module_name,
      "'."
    )
  }
  
  header_index <- header_index[[1]]
  
  header_line <- str_remove(
    module_lines[[header_index]],
    "^#"
  )
  
  if (header_index < length(module_lines)) {
    data_lines <- module_lines[
      seq.int(
        from = header_index + 1L,
        to = length(module_lines)
      )
    ]
  } else {
    data_lines <- character()
  }
  
  data_lines <- data_lines[
    nzchar(str_trim(data_lines))
  ]
  
  if (length(data_lines) == 0L) {
    stop(
      "No data rows were found in FastQC module '",
      module_name,
      "'."
    )
  }
  
  read.delim(
    text = paste(
      c(header_line, data_lines),
      collapse = "\n"
    ),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    comment.char = "",
    quote = "",
    stringsAsFactors = FALSE
  ) %>%
    as_tibble()
}


# Extract one named value from the Basic Statistics module
extract_fastqc_basic_value <- function(
    basic_statistics,
    requested_measure
) {
  
  measure_column <- find_table_column(
    basic_statistics,
    "Measure"
  )
  
  value_column <- find_table_column(
    basic_statistics,
    "Value"
  )
  
  measure_names <- str_to_lower(
    str_trim(
      as.character(
        basic_statistics[[measure_column]]
      )
    )
  )
  
  requested_measure_lower <- str_to_lower(
    requested_measure
  )
  
  selected_row <- which(
    measure_names == requested_measure_lower
  )
  
  if (length(selected_row) == 0L) {
    stop(
      "Measure '",
      requested_measure,
      "' was not found in the FastQC Basic Statistics module."
    )
  }
  
  as.character(
    basic_statistics[[value_column]][selected_row[[1]]]
  )
}


# Prepare one FastQC distribution table
prepare_fastqc_distribution <- function(
    module_table,
    bin_column_name,
    sample_name,
    source_report
) {
  
  bin_column <- find_table_column(
    module_table,
    bin_column_name
  )
  
  count_column <- find_table_column(
    module_table,
    "Count"
  )
  
  tibble(
    Sample = sample_name,
    Source_Report = source_report,
    Bin = as.character(
      module_table[[bin_column]]
    ),
    Value = parse_bin_midpoint(
      module_table[[bin_column]]
    ),
    Count = parse_plain_number(
      module_table[[count_column]]
    )
  ) %>%
    filter(
      is.finite(Value),
      is.finite(Count),
      Count > 0
    )
}


############################################################
# 6) Parse one Illumina FastQC report
############################################################

parse_fastqc_file <- function(
    file,
    data_type
) {
  
  file_name <- basename(file)
  
  # Extract sample name and read direction from the filename
  if (data_type == "Raw") {
    
    sample_name <- str_remove(
      file_name,
      "_raw\\.[12]_fastqc\\.txt$"
    )
    
    read_number <- str_match(
      file_name,
      "_raw\\.([12])_fastqc\\.txt$"
    )[, 2]
    
    if (is.na(read_number)) {
      stop(
        "Could not determine R1 or R2 from raw FastQC filename: ",
        file_name
      )
    }
    
    read_set <- paste0("R", read_number)
    
  } else if (data_type == "QualityFiltered") {
    
    sample_name <- str_remove(
      file_name,
      "\\.qf_fastqc\\.txt$"
    )
    
    read_set <- "Interleaved"
    
  } else {
    stop(
      "Unknown Illumina data type: ",
      data_type
    )
  }
  
  # Load the complete FastQC report
  report_lines <- readLines(
    file,
    warn = FALSE
  )
  
  # Parse the required FastQC modules
  basic_statistics <- read_fastqc_module(
    report_lines,
    "Basic Statistics"
  )
  
  length_module <- read_fastqc_module(
    report_lines,
    "Sequence Length Distribution"
  )
  
  quality_module <- read_fastqc_module(
    report_lines,
    "Per sequence quality scores"
  )
  
  # Prepare the read-length distribution
  length_distribution <- prepare_fastqc_distribution(
    module_table = length_module,
    bin_column_name = "Length",
    sample_name = sample_name,
    source_report = file_name
  )
  
  # Prepare the per-sequence mean-quality distribution
  quality_distribution <- prepare_fastqc_distribution(
    module_table = quality_module,
    bin_column_name = "Quality",
    sample_name = sample_name,
    source_report = file_name
  )
  
  # Extract total number of sequences
  number_of_reads <- extract_fastqc_basic_value(
    basic_statistics,
    "Total Sequences"
  ) %>%
    parse_plain_number()
  
  # Extract total sequence yield.
  #
  # FastQC may report this as, for example, "19.9 Gbp".
  total_bases <- extract_fastqc_basic_value(
    basic_statistics,
    "Total Bases"
  ) %>%
    parse_scaled_bases()
  
  # Calculate mean and median read length from the
  # Sequence Length Distribution module.
  mean_read_length <- weighted_mean_safe(
    length_distribution$Value,
    length_distribution$Count
  )
  
  median_read_length <- weighted_median(
    length_distribution$Value,
    length_distribution$Count
  )
  
  # Calculate mean and median read quality from the
  # Per sequence quality scores module.
  #
  # This is more appropriate than averaging the plotted
  # per-position FastQC means because it represents the
  # distribution of mean quality scores across reads.
  mean_read_quality <- weighted_mean_safe(
    quality_distribution$Value,
    quality_distribution$Count
  )
  
  median_read_quality <- weighted_median(
    quality_distribution$Value,
    quality_distribution$Count
  )
  
  # Store the report-level summary
  report_summary <- tibble(
    Sample = sample_name,
    Technology = "Illumina",
    Data_Type = data_type,
    Read_Set = read_set,
    Mean_Read_Length = mean_read_length,
    Mean_Read_Quality = mean_read_quality,
    Median_Read_Length = median_read_length,
    Median_Read_Quality = median_read_quality,
    Number_Of_Reads = number_of_reads,
    Total_Bases = total_bases,
    Source = "FastQC",
    Source_Report = file_name
  )
  
  # Return all parsed data tables.
  #
  # This preserves the full information needed to pool R1
  # and R2 correctly instead of averaging report summaries.
  list(
    Summary = report_summary,
    Length_Distribution = length_distribution,
    Quality_Distribution = quality_distribution,
    Basic_Statistics = basic_statistics
  )
}


############################################################
# 7) Parse one Nanopore NanoStat report
############################################################

# Extract one labeled value from a NanoStat report
extract_nanostat_value <- function(
    lines,
    label,
    file_name
) {
  
  trimmed_lines <- str_trim(lines)
  
  selected_line <- which(
    startsWith(
      trimmed_lines,
      paste0(label, ":")
    )
  )
  
  if (length(selected_line) == 0L) {
    stop(
      "NanoStat field '",
      label,
      "' was not found in ",
      file_name
    )
  }
  
  value_text <- sub(
    "^[^:]+:",
    "",
    trimmed_lines[selected_line[[1]]]
  )
  
  parse_plain_number(value_text)
}


parse_nanostat_file <- function(file) {
  
  file_name <- basename(file)
  
  # Determine processing status from the filename
  data_type <- case_when(
    str_detect(
      file_name,
      "\\.QualityReport\\.raw\\.txt$"
    ) ~ "Raw",
    
    str_detect(
      file_name,
      "\\.QualityReport\\.qf\\.txt$"
    ) ~ "QualityFiltered",
    
    TRUE ~ NA_character_
  )
  
  if (is.na(data_type)) {
    stop(
      "Could not determine NanoStat data type from filename: ",
      file_name
    )
  }
  
  # Remove the report suffix to recover the sample name
  sample_name <- str_remove(
    file_name,
    "\\.QualityReport\\.(raw|qf)\\.txt$"
  )
  
  report_lines <- readLines(
    file,
    warn = FALSE
  )
  
  report_summary <- tibble(
    Sample = sample_name,
    Technology = "Nanopore",
    Data_Type = data_type,
    Read_Set = "LongReads",
    
    Mean_Read_Length = extract_nanostat_value(
      report_lines,
      "Mean read length",
      file_name
    ),
    
    Mean_Read_Quality = extract_nanostat_value(
      report_lines,
      "Mean read quality",
      file_name
    ),
    
    Median_Read_Length = extract_nanostat_value(
      report_lines,
      "Median read length",
      file_name
    ),
    
    Median_Read_Quality = extract_nanostat_value(
      report_lines,
      "Median read quality",
      file_name
    ),
    
    Number_Of_Reads = extract_nanostat_value(
      report_lines,
      "Number of reads",
      file_name
    ),
    
    Total_Bases = extract_nanostat_value(
      report_lines,
      "Total bases",
      file_name
    ),
    
    Source = "NanoStat",
    Source_Report = file_name
  )
  
  list(
    Summary = report_summary
  )
}


############################################################
# 8) Helper for creating named report lists
############################################################

create_named_report_list <- function(
    files,
    parser,
    ...
) {
  
  if (length(files) == 0L) {
    stop("No input files were found.")
  }
  
  files <- sort(files)
  
  report_list <- lapply(
    files,
    FUN = function(file) {
      parser(file, ...)
    }
  )
  
  report_names <- basename(files) %>%
    str_remove("\\.txt$") %>%
    make.unique()
  
  names(report_list) <- report_names
  
  report_list
}


############################################################
# 9) Find all input reports
############################################################

illumina_raw_files <- list.files(
  path = illumina_raw_dir,
  pattern = "_raw\\.[12]_fastqc\\.txt$",
  full.names = TRUE
)

illumina_filtered_files <- list.files(
  path = illumina_filtered_dir,
  pattern = "\\.qf_fastqc\\.txt$",
  full.names = TRUE
)

nanopore_raw_files <- list.files(
  path = nanopore_dir,
  pattern = "\\.QualityReport\\.raw\\.txt$",
  full.names = TRUE
)

nanopore_filtered_files <- list.files(
  path = nanopore_dir,
  pattern = "\\.QualityReport\\.qf\\.txt$",
  full.names = TRUE
)

cat(
  "Files found:\n",
  "  Illumina raw FastQC:       ",
  length(illumina_raw_files),
  "\n",
  "  Illumina filtered FastQC:  ",
  length(illumina_filtered_files),
  "\n",
  "  Nanopore raw NanoStat:     ",
  length(nanopore_raw_files),
  "\n",
  "  Nanopore filtered NanoStat:",
  length(nanopore_filtered_files),
  "\n\n",
  sep = ""
)


############################################################
# 10) Parse every report individually
############################################################

# Each list is named using the corresponding report filename.
#
# For example:
#
# Illumina_Raw_FastQC_Reports[
#   ["IZE01APR2102_raw.1_fastqc"]
# ]$Summary
#
# returns the parsed summary for that individual report.

Illumina_Raw_FastQC_Reports <- create_named_report_list(
  files = illumina_raw_files,
  parser = parse_fastqc_file,
  data_type = "Raw"
)

Illumina_Filtered_FastQC_Reports <- create_named_report_list(
  files = illumina_filtered_files,
  parser = parse_fastqc_file,
  data_type = "QualityFiltered"
)

Nanopore_Raw_NanoStat_Reports <- create_named_report_list(
  files = nanopore_raw_files,
  parser = parse_nanostat_file
)

Nanopore_Filtered_NanoStat_Reports <- create_named_report_list(
  files = nanopore_filtered_files,
  parser = parse_nanostat_file
)


############################################################
# 11) Create named lists containing only summary tables
############################################################

Illumina_Raw_FastQC_Tables <- lapply(
  Illumina_Raw_FastQC_Reports,
  function(report) report$Summary
)

Illumina_Filtered_FastQC_Tables <- lapply(
  Illumina_Filtered_FastQC_Reports,
  function(report) report$Summary
)

Nanopore_Raw_NanoStat_Tables <- lapply(
  Nanopore_Raw_NanoStat_Reports,
  function(report) report$Summary
)

Nanopore_Filtered_NanoStat_Tables <- lapply(
  Nanopore_Filtered_NanoStat_Reports,
  function(report) report$Summary
)


############################################################
# 12) Check the Illumina report structure
############################################################

illumina_raw_per_file <- bind_rows(
  Illumina_Raw_FastQC_Tables
)

illumina_filtered_per_file <- bind_rows(
  Illumina_Filtered_FastQC_Tables
)

# Every raw Illumina sample should have exactly:
#   one R1 report
#   one R2 report
raw_pair_check <- illumina_raw_per_file %>%
  group_by(Sample) %>%
  summarise(
    Number_Of_Reports = n(),
    Has_R1 = any(Read_Set == "R1"),
    Has_R2 = any(Read_Set == "R2"),
    .groups = "drop"
  )

incorrect_raw_pairs <- raw_pair_check %>%
  filter(
    Number_Of_Reports != 2L |
      !Has_R1 |
      !Has_R2
  )

if (nrow(incorrect_raw_pairs) > 0L) {
  stop(
    "The following raw Illumina samples do not have exactly one R1 ",
    "and one R2 report:\n",
    paste(
      incorrect_raw_pairs$Sample,
      collapse = ", "
    )
  )
}


# Every quality-filtered Illumina sample should have one
# FastQC report representing the interleaved file.
filtered_file_check <- illumina_filtered_per_file %>%
  count(
    Sample,
    name = "Number_Of_Reports"
  ) %>%
  filter(Number_Of_Reports != 1L)

if (nrow(filtered_file_check) > 0L) {
  stop(
    "The following filtered Illumina samples do not have exactly ",
    "one interleaved FastQC report:\n",
    paste(
      filtered_file_check$Sample,
      collapse = ", "
    )
  )
}


############################################################
# 13) Combine Illumina reports at sample level
############################################################

# Raw R1 and R2 are pooled using the underlying FastQC
# distributions. Therefore:
#
# - read counts and total bases are summed;
# - mean read length is weighted by the number of reads;
# - median read length is calculated from the pooled distribution;
# - mean read quality is weighted by the number of reads;
# - median read quality is calculated from the pooled distribution.
#
# The same function also works for filtered reports, for which
# each sample contains one interleaved FastQC report.

combine_fastqc_reports <- function(report_list) {
  
  per_file_summary <- bind_rows(
    lapply(
      report_list,
      function(report) report$Summary
    )
  )
  
  pooled_length_distribution <- bind_rows(
    lapply(
      report_list,
      function(report) report$Length_Distribution
    )
  ) %>%
    group_by(
      Sample,
      Value
    ) %>%
    summarise(
      Count = sum(Count),
      .groups = "drop"
    )
  
  pooled_quality_distribution <- bind_rows(
    lapply(
      report_list,
      function(report) report$Quality_Distribution
    )
  ) %>%
    group_by(
      Sample,
      Value
    ) %>%
    summarise(
      Count = sum(Count),
      .groups = "drop"
    )
  
  base_summary <- per_file_summary %>%
    group_by(
      Sample,
      Technology,
      Data_Type
    ) %>%
    summarise(
      Number_Of_Reads = sum(Number_Of_Reads),
      Total_Bases = sum(Total_Bases),
      Source = first(Source),
      
      Number_Of_Source_Reports = n_distinct(
        Source_Report
      ),
      
      Source_Reports = paste(
        sort(unique(Source_Report)),
        collapse = "; "
      ),
      
      Read_Layout = if_else(
        first(Data_Type) == "Raw",
        "Paired-end (R1 + R2)",
        "Interleaved"
      ),
      
      .groups = "drop"
    )
  
  length_summary <- pooled_length_distribution %>%
    group_by(Sample) %>%
    summarise(
      Mean_Read_Length = weighted_mean_safe(
        Value,
        Count
      ),
      
      Median_Read_Length = weighted_median(
        Value,
        Count
      ),
      
      .groups = "drop"
    )
  
  quality_summary <- pooled_quality_distribution %>%
    group_by(Sample) %>%
    summarise(
      Mean_Read_Quality = weighted_mean_safe(
        Value,
        Count
      ),
      
      Median_Read_Quality = weighted_median(
        Value,
        Count
      ),
      
      .groups = "drop"
    )
  
  base_summary %>%
    left_join(
      length_summary,
      by = "Sample"
    ) %>%
    left_join(
      quality_summary,
      by = "Sample"
    ) %>%
    select(
      Sample,
      Technology,
      Data_Type,
      Read_Layout,
      Mean_Read_Length,
      Mean_Read_Quality,
      Median_Read_Length,
      Median_Read_Quality,
      Number_Of_Reads,
      Total_Bases,
      Source,
      Number_Of_Source_Reports,
      Source_Reports
    )
}


Illumina_Raw_Stats <- combine_fastqc_reports(
  Illumina_Raw_FastQC_Reports
)

Illumina_Filtered_Stats <- combine_fastqc_reports(
  Illumina_Filtered_FastQC_Reports
)


############################################################
# 14) Combine the NanoStat summaries
############################################################

Nanopore_Raw_Stats <- bind_rows(
  Nanopore_Raw_NanoStat_Tables
) %>%
  mutate(
    Read_Layout = "Long reads",
    Number_Of_Source_Reports = 1L,
    Source_Reports = Source_Report
  ) %>%
  select(
    -Read_Set,
    -Source_Report
  )

Nanopore_Filtered_Stats <- bind_rows(
  Nanopore_Filtered_NanoStat_Tables
) %>%
  mutate(
    Read_Layout = "Long reads",
    Number_Of_Source_Reports = 1L,
    Source_Reports = Source_Report
  ) %>%
  select(
    -Read_Set,
    -Source_Report
  )


############################################################
# 15) Extract the sampling date
############################################################

extract_sampling_date <- function(sample_name) {
  
  # Example:
  #
  # IZE01APR2102
  #    01APR21 = 1 April 2021
  #          02 = sampling depth
  #
  # Only DDMMMYY is therefore extracted from the sample name.
  
  date_code <- str_extract(
    sample_name,
    "\\d{2}[A-Z]{3}\\d{2}"
  )
  
  dmy(
    date_code,
    locale = "C",
    quiet = TRUE
  )
}


############################################################
# 16) Create the final unified table
############################################################

QualityFilteringStats <- bind_rows(
  Illumina_Raw_Stats,
  Illumina_Filtered_Stats,
  Nanopore_Raw_Stats,
  Nanopore_Filtered_Stats
) %>%
  mutate(
    Sampling_Date = extract_sampling_date(Sample),
    
    Technology = factor(
      Technology,
      levels = c(
        "Illumina",
        "Nanopore"
      )
    ),
    
    Data_Type = factor(
      Data_Type,
      levels = c(
        "Raw",
        "QualityFiltered"
      )
    )
  ) %>%
  arrange(
    Sampling_Date,
    Technology,
    Data_Type
  ) %>%
  select(
    Sample,
    Sampling_Date,
    Technology,
    Data_Type,
    Read_Layout,
    Mean_Read_Length,
    Mean_Read_Quality,
    Median_Read_Length,
    Median_Read_Quality,
    Number_Of_Reads,
    Total_Bases,
    Source,
    Number_Of_Source_Reports,
    Source_Reports
  )


############################################################
# 17) Final validation
############################################################

# Ensure that every date was recovered successfully
if (any(is.na(QualityFilteringStats$Sampling_Date))) {
  
  problematic_samples <- QualityFilteringStats %>%
    filter(is.na(Sampling_Date)) %>%
    distinct(Sample) %>%
    pull(Sample)
  
  stop(
    "The sampling date could not be extracted for:\n",
    paste(
      problematic_samples,
      collapse = ", "
    )
  )
}


# Ensure that each sample, technology and processing status
# occurs only once in the final table
duplicate_rows <- QualityFilteringStats %>%
  count(
    Technology,
    Data_Type,
    Sample,
    name = "Number_Of_Rows"
  ) %>%
  filter(Number_Of_Rows != 1L)

if (nrow(duplicate_rows) > 0L) {
  stop(
    "Duplicate sample/technology/data-type combinations ",
    "were found in the final table."
  )
}


# Ensure that all principal sequencing metrics were recovered
required_metric_columns <- c(
  "Mean_Read_Length",
  "Mean_Read_Quality",
  "Median_Read_Length",
  "Median_Read_Quality",
  "Number_Of_Reads",
  "Total_Bases"
)

if (
  any(
    !complete.cases(
      QualityFilteringStats[
        required_metric_columns
      ]
    )
  )
) {
  stop(
    "At least one sequencing metric is missing. ",
    "Inspect the relevant source report before exporting."
  )
}


############################################################
# 18) Print a report overview
############################################################

report_overview <- QualityFilteringStats %>%
  count(
    Technology,
    Data_Type,
    name = "Number_Of_Samples"
  )

print(
  report_overview,
  n = Inf
)

cat(
  "\nFinal number of rows: ",
  nrow(QualityFilteringStats),
  "\n",
  sep = ""
)


############################################################
# 19) Export the final CSV file
############################################################

write_csv(
  QualityFilteringStats,
  file = output_file,
  na = "NA"
)

cat(
  "\nQuality-filtering statistics saved to:\n",
  output_file,
  "\n",
  sep = ""
)


############################################################
# 20) Show the completed table in RStudio
############################################################

QualityFilteringStats
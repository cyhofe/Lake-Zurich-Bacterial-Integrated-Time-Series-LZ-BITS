rm(list = ls())
# 1: N2 for each month of 2021 panel plot (March/April/May/June)
    ## Up to 50 m only
    ## Add T; O2; PO4; ChlA profiles
    ## Plot based on 2-3 daily samples
    ## Plot based on median value + inter-4-le range
# 2: N2 dynamics throughout the year for 2-8 m depth
    ## Plot based on data points from 2021, incorporating 2-3 daily samples to the yearly record
    ## Plot each data point from 2021 without marking the dates as dots on X axis, just a line
    ## Each value for 2021 is: N2 calculated for each m within 2-8 meters and then median of it is plotted
    ## For the background: N2 calculated for each m within 2-8 meters for all existing data points for years 2012-2021 and then monthly median with inter-4-le range is plotted, make a smooth line, connecting the months
    ## For the decadal max and min monthly levels use dashed lines


##########################################################################################################################################################################################################################
# 1
################################################################################
# 2021 monthly median profiles (March-June)
# T, O2, PO4 (if present), ChlA, N2
# Line/points = monthly median
# Horizontal bars = interquartile range (Q25-Q75)
# Depth range = 0-50 m
################################################################################

# ------------------------------- packages -------------------------------------
pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "lubridate", "purrr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

invisible(lapply(pkgs, library, character.only = TRUE))

# ------------------------------- paths ----------------------------------------
base_path <- "/SynologyDrive/PhD/Projects/SpringBloom2021/DataAnalysis_SpringBloom2021/Supplementary/LakeProfiles/" 

input_file  <- file.path(base_path, "../../../Data/Lake_Data/Raw_Data/Data_2021.xlsx")
output_csv  <- file.path(base_path, "Profiles_monthly_median_IQR_2021_MarJun.csv")
output_plot <- file.path(base_path, "Profiles_monthly_median_IQR_2021_MarJun_panel.png")

# ------------------------------ helpers ---------------------------------------
to_num <- function(x) {
  if (is.numeric(x)) return(x)
  x <- trimws(as.character(x))
  x[x == ""] <- NA
  as.numeric(gsub(",", ".", x, fixed = TRUE))
}

parse_excel_date <- function(x) {
  if (inherits(x, "Date"))   return(as.Date(x))
  if (inherits(x, "POSIXt")) return(as.Date(x))
  
  x_chr <- as.character(x)
  out   <- as.Date(rep(NA_character_, length(x_chr)))
  
  # text dates like 2021-02-22
  d_txt <- suppressWarnings(lubridate::ymd(x_chr))
  
  # Excel serial dates if needed
  d_num <- suppressWarnings(as.Date(as.numeric(x_chr), origin = "1899-12-30"))
  
  idx_txt <- !is.na(d_txt)
  out[idx_txt] <- as.Date(d_txt[idx_txt])
  
  idx_num <- is.na(out) & !is.na(d_num)
  out[idx_num] <- d_num[idx_num]
  
  out
}

q_safe <- function(x, p) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  unname(quantile(x, probs = p, type = 7))
}

find_first_col <- function(nms, candidates) {
  hit <- nms[nms %in% candidates]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# Freshwater density from temperature (kg m-3)
rho_fresh <- function(temp_c) {
  999.842594 +
    6.793952e-2 * temp_c -
    9.095290e-3 * temp_c^2 +
    1.001685e-4 * temp_c^3 -
    1.120083e-6 * temp_c^4 +
    6.536332e-9 * temp_c^5
}

# N2 for one profile
# depth is positive downward, so:
# N2 = g / rho * d(rho)/dz
calc_n2_one_profile <- function(df_one_profile) {
  g <- 9.81
  
  x <- df_one_profile %>%
    dplyr::filter(!is.na(depth), !is.na(temperature)) %>%
    dplyr::arrange(depth) %>%
    dplyr::distinct(depth, .keep_all = TRUE)
  
  if (nrow(x) < 2) {
    return(tibble(
      date      = as.Date(character()),
      month_lab = factor(character(), levels = c("March", "April", "May", "June")),
      depth     = numeric(),
      parameter = character(),
      value     = numeric()
    ))
  }
  
  rho  <- rho_fresh(x$temperature)
  dz   <- diff(x$depth)
  drho <- diff(rho)
  
  ok <- !is.na(dz) & dz > 0 & !is.na(drho)
  
  if (!any(ok)) {
    return(tibble(
      date      = as.Date(character()),
      month_lab = factor(character(), levels = c("March", "April", "May", "June")),
      depth     = numeric(),
      parameter = character(),
      value     = numeric()
    ))
  }
  
  n2 <- (g / ((rho[-1] + rho[-length(rho)]) / 2)) * (drho / dz)
  
  tibble(
    date      = x$date[1],
    month_lab = x$month_lab[1],
    depth     = (x$depth[-1] + x$depth[-nrow(x)]) / 2,
    parameter = "N2",
    value     = n2
  ) %>%
    filter(!is.na(value))
}

# ------------------------------- read data ------------------------------------
dat <- read_excel(input_file)

# clean names
names(dat) <- tolower(gsub("[^[:alnum:]]+", "_", names(dat)))
names(dat) <- gsub("^_|_$", "", names(dat))

# detect main columns
date_col <- find_first_col(names(dat), c("date"))
depth_col <- find_first_col(names(dat), c("depth"))
temp_col <- find_first_col(names(dat), c("temperature", "temp"))
o2_col   <- find_first_col(names(dat), c("oxygen", "o2"))
chl_col  <- find_first_col(names(dat), c("chlorophyll", "chla", "chl_a", "chlorophyll_a"))
po4_col  <- find_first_col(names(dat), c("po4", "phosphate", "srp"))

required_cols <- c(date_col, depth_col, temp_col, o2_col, chl_col)
if (any(is.na(required_cols))) {
  stop(
    "Missing one or more required columns.\n",
    "I need at least: date, depth, temperature, oxygen, chlorophyll.\n",
    "Detected columns were:\n",
    paste(names(dat), collapse = ", ")
  )
}

# standardize names used later
dat <- dat %>%
  mutate(
    date        = parse_excel_date(.data[[date_col]]),
    depth       = to_num(.data[[depth_col]]),
    temperature = to_num(.data[[temp_col]]),
    oxygen      = to_num(.data[[o2_col]]),
    chlorophyll = to_num(.data[[chl_col]])
  )

if (!is.na(po4_col)) {
  dat$po4 <- to_num(dat[[po4_col]])
}

# keep only March-June 2021 and top 50 m
dat_2021 <- dat %>%
  mutate(
    year      = lubridate::year(date),
    month     = lubridate::month(date),
    month_lab = factor(
      month.name[month],
      levels = c("March", "April", "May", "June")
    )
  ) %>%
  filter(
    year == 2021,
    month %in% 3:6,
    !is.na(depth),
    depth >= 0,
    depth <= 50
  )

# ------------------------ build long table: measured vars ----------------------
profiles_list <- list(
  T = dat_2021 %>%
    transmute(date, month_lab, depth, parameter = "T", value = temperature),
  
  O2 = dat_2021 %>%
    transmute(date, month_lab, depth, parameter = "O2", value = oxygen),
  
  ChlA = dat_2021 %>%
    transmute(date, month_lab, depth, parameter = "ChlA", value = chlorophyll)
)

if ("po4" %in% names(dat_2021)) {
  profiles_list$PO4 <- dat_2021 %>%
    transmute(date, month_lab, depth, parameter = "PO4", value = po4)
}

# ----------------------------- N2 profiles ------------------------------------
n2_profiles <- dat_2021 %>%
  group_by(date, month_lab) %>%
  group_split() %>%
  purrr::map_dfr(calc_n2_one_profile)

# ----------------------- combine all parameters --------------------------------
profiles_all <- bind_rows(
  bind_rows(profiles_list),
  n2_profiles
) %>%
  filter(!is.na(value))

param_order <- c("T", "O2", "PO4", "ChlA", "N2")
param_order <- param_order[param_order %in% unique(profiles_all$parameter)]

profiles_all <- profiles_all %>%
  mutate(
    parameter = factor(parameter, levels = param_order)
  )

# ------------------------ monthly median + IQR ---------------------------------
monthly_summary <- profiles_all %>%
  group_by(parameter, month_lab, depth) %>%
  summarise(
    median_value = median(value, na.rm = TRUE),
    q25_value    = q_safe(value, 0.25),
    q75_value    = q_safe(value, 0.75),
    n_profiles   = sum(!is.na(value)),
    .groups = "drop"
  ) %>%
  arrange(parameter, month_lab, depth) %>%
  filter(!is.na(median_value))

write.csv(monthly_summary, output_csv, row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------------- plot ------------------------------------------
month_cols <- c(
  "March" = "#1b9e77",
  "April" = "#7570b3",
  "May"   = "#d95f02",
  "June"  = "#e7298a"
)

facet_labels <- c(
  "T"    = "T",
  "O2"   = "O2",
  "PO4"  = "PO4",
  "ChlA" = "ChlA",
  "N2"   = "N2"
)

p <- ggplot(
  monthly_summary,
  aes(x = median_value, y = depth, colour = month_lab, group = month_lab)
) +
  geom_segment(
    aes(x = q25_value, xend = q75_value, y = depth, yend = depth),
    linewidth = 0.35,
    alpha = 0.50,
    show.legend = FALSE
  ) +
  geom_path(linewidth = 0.9, na.rm = TRUE) +
  geom_point(size = 1.4, na.rm = TRUE) +
  geom_hline(
    yintercept = c(2, 8),
    linetype = "dashed",
    colour = "grey40"
  ) +
  scale_y_reverse(
    limits = c(50, 0),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_colour_manual(
    values = month_cols,
    drop = FALSE,
    name = "Month"
  ) +
  facet_wrap(
    ~ parameter,
    scales = "free_x",
    ncol = 2,
    labeller = as_labeller(facet_labels)
  ) +
  labs(
    x = NULL,
    y = "Depth (m)",
    title = "2021 monthly median vertical profiles",
    subtitle = "March-June; horizontal bars show interquartile range (Q25-Q75)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )

print(p)

ggsave(
  filename = output_plot,
  plot = p,
  width = 11,
  height = 8.5,
  dpi = 300
)


# ============================================================
# Source data export: Supplementary Figure 3
# One row = one date × depth × parameter value
# ============================================================

base_dir <- paste0(
  "C:/SynologyDrive/PhD/Projects/SpringBloom2021/",
  "DataAnalysis_SpringBloom2021"
)

source_dir <- file.path(
  base_dir,
  "Source_Data_Supplement"
)

dir.create(
  source_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

suppl_figure_3 <- profiles_all |>
  dplyr::transmute(
    Date = as.Date(date),
    Month = as.character(month_lab),
    Depth_m = as.numeric(depth),
    Parameter = as.character(parameter),
    Observation_Type = dplyr::if_else(
      as.character(parameter) == "N2",
      "Calculated from temperature profile",
      "Measured"
    ),
    Parameter_Value = as.numeric(value)
  ) |>
  dplyr::arrange(
    Parameter,
    Month,
    Date,
    Depth_m
  )

readr::write_excel_csv(
  suppl_figure_3,
  file.path(source_dir, "Suppl_Figure_3.csv"),
  na = ""
)

profiles_all %>%
  distinct(month_lab, date) %>%
  count(month_lab, name = "n_sampling_dates") %>%
  print(n = Inf)

profiles_all %>%
  distinct(parameter, month_lab, date) %>%
  count(parameter, month_lab, name = "n_sampling_dates") %>%
  print(n = Inf)

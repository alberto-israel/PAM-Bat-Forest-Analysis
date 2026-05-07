# ============================================================
# 07_temperature.R
# Author: Alberto Israel
# Purpose: Extract temperature data from Kaleidoscope summary files
#          and produce hourly and daily temperature tables
# Input:   data/raw_pam_outputs/S*.txt (Wildlife Acoustics / Kaleidoscope summary)
# Output:  outputs/temperature_site.csv
#          outputs/temperature_daily.csv
#          outputs/temperature_hourly.csv
# Notes:
#   - Requires summary .txt files in data/raw_pam_outputs/
#     (naming convention: S1.txt, S2.txt, ... SN.txt)
#   - Nocturnal sessions defined as 18:00–08:00 (consistent with 00b)
#   - Joinable with the site × species matrices created in 00b
#   - Run after 00b_data_preparation.R and before analysis scripts
# Last update: 2026-05-06
# ============================================================


# ------------------------------------------------------------
# 0) Setup project environment
# ------------------------------------------------------------
source(here::here("scripts", "00a_setup.R"))

data_dir   <- file.path(here::here(), "data")
output_dir <- file.path(here::here(), "outputs")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# ============================================================
# 1) Import summary files
# ============================================================

# ------------------------------------------------------------
# 1.1) Locate summary files
# ------------------------------------------------------------

# Summary files are expected directly inside raw_pam_outputs/
# Expected naming: S1.txt, S2.txt, ..., SN.txt
# If your files follow a different convention, adjust the pattern below

raw_data_dir  <- file.path(data_dir, "raw_pam_outputs")
summary_files <- list.files(raw_data_dir, pattern = "^S[0-9]+\\.txt$", full.names = TRUE)

cat("Found", length(summary_files), "summary file(s)\n")

if (length(summary_files) == 0) {
  stop(
    "No summary files found in ", raw_data_dir,
    "\nExpected naming: S1.txt, S2.txt, ..., SN.txt"
  )
}

# ------------------------------------------------------------
# 1.2) Import and combine all files
# ------------------------------------------------------------

# Kaleidoscope summary format (comma-separated):
# DATE, TIME, LAT, NS, LON, EW, POWER(V), TEMP(C), #FSFILES, #ZCFILES, #SCRUBBED
# Resolution: 1 row per minute

temp_raw <- map_dfr(summary_files, function(f) {
  site_id <- str_extract(basename(f), "[0-9]+") |> as.integer()
  read_csv(f, show_col_types = FALSE) |>
    mutate(SiteID = site_id)
})

cat("Total records imported:", nrow(temp_raw), "\n")
cat("Sites detected:", length(unique(temp_raw$SiteID)), "\n\n")


# ============================================================
# 2) Parse and clean
# ============================================================

temp_clean <- temp_raw |>
  rename(
    date_raw = DATE,
    time_raw = TIME,
    lat      = LAT,
    lon      = LON,
    temp_c   = `TEMP(C)`
  ) |>
  mutate(
    date = as.Date(date_raw, format = "%Y-%b-%d"),  # e.g. 2024-Jun-13
    hour = as.integer(substr(time_raw, 1, 2))       # extract HH from HH:MM:SS
  ) |>
  filter(!is.na(temp_c), !is.na(date), !is.na(hour))

# Check and report missing values removed
n_removed <- nrow(temp_raw) - nrow(temp_clean)
if (n_removed > 0) {
  cat("Records removed (missing date/temp):", n_removed, "\n")
}

# ------------------------------------------------------------
# 2.1) Assign nocturnal sessions
# ------------------------------------------------------------

# Nocturnal session definition (consistent with 00b_data_preparation.R):
#   - Evening: HOUR >= 18 → session date = calendar date
#   - Morning: HOUR < 8   → session date = calendar date - 1 day
#   - Daytime records (08:00–18:00) are retained but flagged

temp_clean <- temp_clean |>
  mutate(
    period = case_when(
      hour < 8  ~ "morning",
      hour >= 18 ~ "evening",
      TRUE       ~ "daytime"
    ),
    session_date = if_else(period == "morning", date - days(1), date)
  )

cat("Session boundaries: 18:00 (evening) to 08:00 next day (morning)\n")
cat("Daytime records (08:00–18:00) retained but excluded from nocturnal aggregations\n\n")


# ============================================================
# 3) Site-level temperature table
# ============================================================
# Mean, min, max, and SD of nocturnal temperature per site
# Aggregated across all nocturnal records (18:00–08:00)

temperature_site <- temp_clean |>
  filter(period != "daytime") |>
  group_by(SiteID) |>
  summarise(
    mean_temp_nocturnal = round(mean(temp_c, na.rm = TRUE), 2),
    min_temp_nocturnal  = round(min(temp_c,  na.rm = TRUE), 2),
    max_temp_nocturnal  = round(max(temp_c,  na.rm = TRUE), 2),
    sd_temp_nocturnal   = round(sd(temp_c,   na.rm = TRUE), 2),
    n_nights            = n_distinct(session_date),
    .groups = "drop"
  )

write_csv(temperature_site, file.path(output_dir, "temperature_site.csv"))
cat("Site-level temperature table saved to outputs/temperature_site.csv\n")
cat("Dimensions:", nrow(temperature_site), "rows ×", ncol(temperature_site), "columns\n\n")


# ============================================================
# 4) Daily temperature table
# ============================================================

# Mean, min, max temperature per site × nocturnal session
# Aggregated across all nocturnal hours (18:00–08:00)

temperature_daily <- temp_clean |>
  filter(period != "daytime") |>
  group_by(SiteID, session_date) |>
  summarise(
    mean_temp   = mean(temp_c, na.rm = TRUE),
    min_temp    = min(temp_c,  na.rm = TRUE),
    max_temp    = max(temp_c,  na.rm = TRUE),
    temp_range  = max_temp - min_temp,
    n_records   = n(),
    .groups     = "drop"
  ) |>
  rename(Day = session_date) |>
  arrange(SiteID, Day)

write_csv(temperature_daily, file.path(output_dir, "temperature_daily.csv"))
cat("Daily temperature table saved to outputs/temperature_daily.csv\n")
cat("Dimensions:", nrow(temperature_daily), "rows ×", ncol(temperature_daily), "columns\n\n")


# ============================================================
# 5) Hourly temperature table
# ============================================================

# Mean temperature per site × nocturnal session × hour
# Nocturnal only (evening + morning); daytime excluded

temperature_hourly <- temp_clean |>
  filter(period != "daytime") |>
  group_by(SiteID, session_date, hour) |>
  summarise(
    mean_temp   = mean(temp_c, na.rm = TRUE),
    min_temp    = min(temp_c,  na.rm = TRUE),
    max_temp    = max(temp_c,  na.rm = TRUE),
    n_records   = n(),
    .groups     = "drop"
  ) |>
  rename(Day = session_date, Hour = hour) |>
  arrange(SiteID, Day, Hour)

write_csv(temperature_hourly, file.path(output_dir, "temperature_hourly.csv"))
cat("Hourly temperature table saved to outputs/temperature_hourly.csv\n")
cat("Dimensions:", nrow(temperature_hourly), "rows ×", ncol(temperature_hourly), "columns\n\n")


# ============================================================
# 6) Quality checks
# ============================================================

cat("--- Temperature Quality Report ---\n")
cat("Date range:", format(min(temp_clean$date)), "to", format(max(temp_clean$date)), "\n")
cat("Temperature range across all sites:",
    round(min(temp_clean$temp_c, na.rm = TRUE), 2), "°C to",
    round(max(temp_clean$temp_c, na.rm = TRUE), 2), "°C\n\n")

# ------------------------------------------------------------
# 6.1) Per-site summary
# ------------------------------------------------------------
site_temp_summary <- temp_clean |>
  group_by(SiteID) |>
  summarise(
    n_nights    = n_distinct(session_date[period != "daytime"]),
    mean_temp   = round(mean(temp_c, na.rm = TRUE), 2),
    min_temp    = round(min(temp_c,  na.rm = TRUE), 2),
    max_temp    = round(max(temp_c,  na.rm = TRUE), 2),
    .groups     = "drop"
  )

print(site_temp_summary)
cat("\n")

# ------------------------------------------------------------
# 6.2) Flag sites with suspiciously low record counts
# ------------------------------------------------------------
low_records <- temperature_daily |>
  filter(n_records < 60) |>
  dplyr::select(SiteID, Day, n_records)

if (nrow(low_records) > 0) {
  cat("WARNING: Nights with fewer than 60 minute-records (possible recording gaps):\n")
  print(low_records)
} else {
  cat("No recording gaps detected.\n")
}
cat("\n")


# ============================================================
# 7) Optional: join temperature with detection matrices
# ============================================================
# Requires outputs from 00b_data_preparation.R to be present in outputs/.

# ------------------------------------------------------------
# 7.1) Join with site x species matrices
# ------------------------------------------------------------
# Default: reads from outputs/ (generated by 00b_data_preparation.R)
# Option B users: either copy your matrix to outputs/ or change the path below

bat_data <- read_csv(file.path(output_dir, "bat_data_site_species.csv"),
                     show_col_types = FALSE) |>
  mutate(SiteID = as.integer(SiteID))

bat_data <- bat_data |>
  left_join(temperature_site, by = "SiteID")

write_csv(bat_data, file.path(output_dir, "bat_data_site_species.csv"))


# ------------------------------------------------------------
# 7.2) Join with daily site x specie matrices
# ------------------------------------------------------------
bat_data_daily <- read_csv(file.path(output_dir, "bat_data_daily.csv"),
                            show_col_types = FALSE) |>
   mutate(SiteID = as.integer(SiteID), Day = as.Date(Day))

bat_data_daily <- bat_data_daily |>
   left_join(temperature_daily, by = c("SiteID", "Day"))

write_csv(bat_data_daily, file.path(output_dir, "bat_data_daily.csv"))


# ------------------------------------------------------------
# 7.3) Join with hourly site x specie matrices
# ------------------------------------------------------------
bat_data_hourly <- read_csv(file.path(output_dir, "bat_data_hourly.csv"),
                            show_col_types = FALSE) |>
  mutate(SiteID = as.integer(SiteID), Day = as.Date(Day), Hour = as.integer(Hour))

bat_data_hourly <- bat_data_hourly |>
  left_join(temperature_hourly, by = c("SiteID", "Day", "Hour"))

write_csv(bat_data_hourly, file.path(output_dir, "bat_data_hourly.csv"))


# ============================================================
# 8) Final output summary
# ============================================================

cat("--- Files created in outputs/ ---\n")
cat("  - temperature_site.csv     (SiteID: mean nocturnal temp across all nights)\n")
cat("  - temperature_daily.csv    (SiteID × Day: mean/min/max/range temp)\n")
cat("  - temperature_hourly.csv   (SiteID × Day × Hour: mean/min/max temp)\n")
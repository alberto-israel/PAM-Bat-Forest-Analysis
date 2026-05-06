# ============================================================
# 00b_data_preparation.R
# Author: Alberto Israel
# Purpose: Transform raw PAM output files into analysis-ready matrices
# Methods: Multi-file import, data validation, temporal aggregation
# Notes:
#   - Supports Wildlife Acoustics / Kaleidoscope output format
#   - Processes multiple sites in batch
#   - Creates multiple aggregation levels (site, day, hour)
#   - Extracts site coordinates from Kaleidoscope summary files (S*.txt)
#   - OPTIONAL: Only needed if starting from raw Kaleidoscope output
#   - Skip if you already have aggregated site × species matrix
#   - Designed to precede 01_data_import.R
# Last update: 2026-05-06
# ============================================================

# ------------------------------------------------------------
# 0) Setup project environment
# ------------------------------------------------------------
source(here::here("scripts", "00a_setup.R"))

# Define directories
data_dir <- file.path(here::here(), "data")
output_dir <- file.path(here::here(), "outputs")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================
# 1) Multi-site Import
# ============================================================

# ------------------------------------------------------------
# 1.1) Define input files
# ------------------------------------------------------------

# NOTE:
# Kaleidoscope creates one output folder per site, each containing id.csv
# Expected structure:
#   data/raw_pam_outputs/
#   ├── Output_Punto_1/
#   │   ├── id.csv
#   │   └── [other files]
#   ├── Output_Punto_2/
#   │   ├── id.csv
#   │   └── [other files]
#   └── Output_Punto_N/
#       ├── id.csv
#       └── [other files]

# switch directory in the script if needed
raw_data_dir <- file.path(data_dir, "raw_pam_outputs")

# Find all output folders (subdirectories)
output_folders <- list.dirs(raw_data_dir, recursive = FALSE, full.names = TRUE)

cat("Found", length(output_folders), "output folders\n")

if (length(output_folders) == 0) {
  stop("No output folders found in ", raw_data_dir,
       "\nPlease place your Kaleidoscope output folders there.")
}

# Locate id.csv in each folder
csv_files <- file.path(output_folders, "id.csv")

# Check which files actually exist
csv_files_exist <- file.exists(csv_files)
csv_files <- csv_files[csv_files_exist]

cat("Found", length(csv_files), "id.csv files to process\n")

if (length(csv_files) == 0) {
  stop("No id.csv files found in output folders.",
       "\nCheck that each folder contains an id.csv file.")
}

# ------------------------------------------------------------
# 1.2) Import and combine all files
# ------------------------------------------------------------

# Function to import single file and extract site ID
import_single_file <- function(filepath) {

  # Read CSV
  data <- read_csv(filepath, show_col_types = FALSE)

  # Extract site ID from parent folder name
  # Example: "...raw_pam_outputs/Output_Punto_14/id.csv" → "14"
  folder_name <- basename(dirname(filepath))

  # Try multiple extraction patterns
  site_id <- str_extract(folder_name, "\\d+")  # Extract first number

  # If extraction fails, try from INDIR column as fallback
  if (is.na(site_id) & "INDIR" %in% colnames(data)) {
    indir_path <- data$INDIR[1]
    site_id <- str_extract(indir_path, "Punto\\s*(\\d+)")[1]
    site_id <- str_extract(site_id, "\\d+")
  }

  # If still NA, use folder name as-is
  if (is.na(site_id)) {
    site_id <- folder_name
    warning("Could not extract numeric site ID from ", folder_name,
            ". Using folder name as Site ID.")
  }

  # Add site ID column
  data$Site <- site_id

  return(data)
}

# Import all files and combine
message("Importing files...")

bat_raw <- map_dfr(csv_files, import_single_file)

cat("Total records imported:", nrow(bat_raw), "\n")
cat("Sites detected:", length(unique(bat_raw$Site)), "\n\n")


# ============================================================
# 2) Data Validation & Quality Control
# ============================================================

# ------------------------------------------------------------
# 2.1) Check for required columns
# ------------------------------------------------------------

required_cols <- c("DATE", "TIME", "HOUR", "MANUAL ID*", "Site")

missing_cols <- setdiff(required_cols, colnames(bat_raw))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("All required columns present\n")

# ------------------------------------------------------------
# 2.2) Remove noise and empty records
# ------------------------------------------------------------

# Count records before filtering
n_before <- nrow(bat_raw)

# Remove noise
bat_clean <- bat_raw %>%
  filter(`MANUAL ID*` != "Noise" | is.na(`MANUAL ID*`))

# Remove empty species IDs
bat_clean <- bat_clean %>%
  filter(!is.na(`MANUAL ID*`), `MANUAL ID*` != "")

n_after <- nrow(bat_clean)

cat("Records removed (Noise/Empty):", n_before - n_after, "\n")
cat("Records retained:", n_after, "\n\n")

# ------------------------------------------------------------
# 2.3) Date/time validation
# ------------------------------------------------------------

# Convert to proper date/time formats
bat_clean <- bat_clean %>%
  mutate(
    DATE = as.Date(DATE),
    TIME = hms(TIME),
    HOUR = as.integer(HOUR)
  )

# Check for invalid dates
invalid_dates <- sum(is.na(bat_clean$DATE))

if (invalid_dates > 0) {
  warning("Found ", invalid_dates, " records with invalid dates. Review manually.")
}

# ------------------------------------------------------------
# 2.4) Species name standardization
# ------------------------------------------------------------

# Standardize species codes: uppercase and trim whitespace
bat_clean <- bat_clean %>%
  mutate(Species = str_trim(toupper(`MANUAL ID*`)))

# Report unique species detected
cat("Unique species detected:", length(unique(bat_clean$Species)), "\n")
cat(paste(sort(unique(bat_clean$Species)), collapse = ", "), "\n\n")


# ============================================================
# 3) Temporal Aggregation - Nocturnal Sessions
# ============================================================

# ------------------------------------------------------------
# 3.1) Define session boundaries
# ------------------------------------------------------------

# NOTE:
# Bat activity is nocturnal. A "session" spans from evening to next morning.
# Example: Session 2024-06-10 includes:
#   - 2024-06-10 from 18:00 to 23:59
#   - 2024-06-11 from 00:00 to 08:00

# Define morning/evening based on hour
bat_clean <- bat_clean %>%
  mutate(
    period = case_when(
      HOUR < 8 ~ "morning",
      HOUR >= 18 ~ "evening",
      TRUE ~ "daytime"
    )
  )

# Assign session date
# Morning records belong to previous day's session
bat_clean <- bat_clean %>%
  mutate(
    session_date = if_else(period == "morning",
                           DATE - days(1),
                           DATE)
  )

cat("Sessions defined from 18:00 to 08:00 next day\n")

# ------------------------------------------------------------
# 3.2) Filter out daytime records (optional)
# ------------------------------------------------------------

# Remove records outside bat activity hours (08:00-18:00)
daytime_records <- sum(bat_clean$period == "daytime")

if (daytime_records > 0) {
  cat("Removing", daytime_records, "daytime records (08:00-18:00)\n")
  bat_clean <- bat_clean %>%
    filter(period != "daytime")
}


# ============================================================
# 4) Create Analysis Matrices
# ============================================================

# ------------------------------------------------------------
# 4.1) Site × Species Matrix (aggregated)
# ------------------------------------------------------------

# Count total calls per species per site
matrix_site_species <- bat_clean %>%
  group_by(Site, Species) %>%
  summarize(N_calls = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = Species,
    values_from = N_calls,
    values_fill = 0
  )

# Add summary columns
matrix_site_species <- matrix_site_species %>%
  mutate(N_calls_total = rowSums(select(., -Site)))

matrix_site_species <- matrix_site_species %>%
  mutate(
    N_species = rowSums(select(., -Site, -N_calls_total) > 0),
    recording_days = NA  # Placeholder - fill manually or from metadata
  )

# Save
write_csv(matrix_site_species,
          file.path(output_dir, "bat_data_site_species.csv"))

cat("  Saved:", file.path(output_dir, "bat_data_site_species.csv"), "\n")
cat("  Dimensions:", nrow(matrix_site_species), "sites ×",
    ncol(matrix_site_species) - 3, "species\n\n")

# ------------------------------------------------------------
# 4.2) Site × Day × Species Matrix
# ------------------------------------------------------------

# Count calls per species per site per day
matrix_site_day_species <- bat_clean %>%
  group_by(Site, session_date, Species) %>%
  summarize(N_calls = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = Species,
    values_from = N_calls,
    values_fill = 0
  ) %>%
  rename(Day = session_date)

# Add summary columns (two-step)
matrix_site_day_species <- matrix_site_day_species %>%
  mutate(N_calls_total = rowSums(select(., -Site, -Day)))

matrix_site_day_species <- matrix_site_day_species %>%
  mutate(N_species = rowSums(select(., -Site, -Day, -N_calls_total) > 0))

# Save
write_csv(matrix_site_day_species,
          file.path(output_dir, "bat_data_daily.csv"))

cat("  Saved:", file.path(output_dir, "bat_data_daily.csv"), "\n")
cat("  Dimensions:", nrow(matrix_site_day_species), "site-days ×",
    ncol(matrix_site_day_species) - 4, "species\n\n")

# ------------------------------------------------------------
# 4.3) Site × Hour × Day × Species Matrix
# ------------------------------------------------------------

# Define hourly bins
matrix_site_hour_species <- bat_clean %>%
  group_by(Site, session_date, HOUR, Species) %>%
  summarize(N_calls = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = Species,
    values_from = N_calls,
    values_fill = 0
  ) %>%
  rename(Day = session_date, Hour = HOUR)

# Add summary columns (two-step)
matrix_site_hour_species <- matrix_site_hour_species %>%
  mutate(N_calls_total = rowSums(select(., -Site, -Day, -Hour)))

matrix_site_hour_species <- matrix_site_hour_species %>%
  mutate(N_species = rowSums(select(., -Site, -Day, -Hour, -N_calls_total) > 0))

# Save
write_csv(matrix_site_hour_species,
          file.path(output_dir, "bat_data_hourly.csv"))

cat("  Saved:", file.path(output_dir, "bat_data_hourly.csv"), "\n")
cat("  Dimensions:", nrow(matrix_site_hour_species), "site-hour-days ×",
    ncol(matrix_site_hour_species) - 5, "species\n\n")


# ============================================================
# 5) Extract Site Coordinates
# ============================================================

# Coordinates (LAT, LON) are not available in id.csv files.
# They are extracted from Kaleidoscope summary files (S*.txt),
# which are expected in data/raw_pam_outputs/ alongside the output folders.
# Expected naming: S1.txt, S2.txt, ..., SN.txt

summary_files <- list.files(raw_data_dir, pattern = "^S[0-9]+\\.txt$",
                             full.names = TRUE)

if (length(summary_files) == 0) {
  warning(
    "No summary files (S*.txt) found in ", raw_data_dir,
    "\nSite coordinates could not be extracted.",
    "\nsite_coordinates.csv will not be created.",
    "\nPlace Kaleidoscope summary files in raw_pam_outputs/ and re-run."
  )
} else {
  cat("Found", length(summary_files), "summary file(s) for coordinate extraction\n")

  coords_raw <- map_dfr(summary_files, function(f) {
    site_id <- str_extract(basename(f), "[0-9]+") |> as.integer()
    read_csv(f, show_col_types = FALSE) |>
      dplyr::select(LAT, LON) |>
      slice(1) |>
      mutate(SiteID = site_id)
  })

  site_coords <- coords_raw |>
    rename(lat = LAT, lon = LON) |>
    dplyr::select(SiteID, lat, lon) |>
    arrange(SiteID)

  write_csv(site_coords, file.path(output_dir, "site_coordinates.csv"))
  cat("Site coordinates saved to outputs/site_coordinates.csv\n")
  print(site_coords)
  cat("\n")
}


# ============================================================
# 6) Environmental Data Integration (if available)
# ============================================================

# If you have habitat classification or environmental data:
# - Load CSV with Site → Habitat mapping
# - Join with matrices

# Example:
# habitat_data <- read_csv(file.path(data_dir, "site_habitat.csv"))
#
# matrix_site_species <- matrix_site_species %>%
#   left_join(habitat_data, by = "Site")

cat("Load habitat/environmental data manually if available.\n\n")


# ============================================================
# 7) Data Quality Report
# ============================================================

# ------------------------------------------------------------
# 7.1) Summary statistics
# ------------------------------------------------------------

cat("--- Summary Statistics ---\n")
cat("Total sites:", length(unique(bat_clean$Site)), "\n")
cat("Total sampling sessions:", length(unique(bat_clean$session_date)), "\n")
cat("Date range:", format(min(bat_clean$DATE)), "to", format(max(bat_clean$DATE)), "\n")
cat("Total bat calls:", nrow(bat_clean), "\n")
cat("Species detected:", length(unique(bat_clean$Species)), "\n\n")

# ------------------------------------------------------------
# 7.2) Per-site summary
# ------------------------------------------------------------

site_summary <- bat_clean %>%
  group_by(Site) %>%
  summarize(
    N_calls = n(),
    N_species = n_distinct(Species),
    N_days = n_distinct(session_date),
    Date_start = min(DATE),
    Date_end = max(DATE)
  )

print(site_summary)

# Save site summary
write_csv(site_summary,
          file.path(output_dir, "site_summary.csv"))

cat("\nSite summary saved to:", file.path(output_dir, "site_summary.csv"), "\n")

# ------------------------------------------------------------
# 7.3) Species frequency
# ------------------------------------------------------------

species_summary <- bat_clean %>%
  group_by(Species) %>%
  summarize(
    N_calls = n(),
    N_sites = n_distinct(Site),
    Frequency = n() / nrow(bat_clean) * 100
  ) %>%
  arrange(desc(N_calls))

print(species_summary)

# Save species summary
write_csv(species_summary,
          file.path(output_dir, "species_frequency.csv"))

cat("\nSpecies summary saved to:", file.path(output_dir, "species_frequency.csv"), "\n")


# ============================================================
# 8) Validation Checks
# ============================================================

# ------------------------------------------------------------
# 8.1) Check for duplicates
# ------------------------------------------------------------

duplicates <- bat_clean %>%
  group_by(Site, DATE, TIME, Species) %>%
  filter(n() > 1) %>%
  nrow()

if (duplicates > 0) {
  warning("Found ", duplicates, " potential duplicate records. Review manually.")
} else {
  cat("No duplicates detected.\n")
}

# ------------------------------------------------------------
# 8.2) Check for suspicious patterns
# ------------------------------------------------------------

# Sites with very low activity
low_activity_sites <- site_summary %>%
  filter(N_calls < 10) %>%
  pull(Site)

if (length(low_activity_sites) > 0) {
  cat("WARNING: Sites with <10 calls:",
      paste(low_activity_sites, collapse = ", "), "\n")
  cat("Review data quality for these sites.\n")
}

# Species with very low occurrence
rare_species <- species_summary %>%
  filter(N_calls < 5) %>%
  pull(Species)

if (length(rare_species) > 0) {
  cat("NOTE: Rare species (<5 calls):",
      paste(rare_species, collapse = ", "), "\n")
}

cat("\n")

# ============================================================
# 9) Final Output Summary
# ============================================================

output_files <- c(
  "bat_data_site_species.csv",
  "bat_data_daily.csv",
  "bat_data_hourly.csv",
  "site_coordinates.csv",
  "site_summary.csv",
  "species_frequency.csv"
)

cat("The following files have been created in outputs/:\n")
for (f in output_files) {
  cat("  -", f, "\n")
}

cat("\n--- Next Steps ---\n")
cat("1. Review site_summary.csv and species_frequency.csv for data quality\n")
cat("2. If you have temperature data (S*.txt), run 07_temperature.R\n")
cat("3. For lunar covariates, run 08_lunar_covariates.R\n")
cat("4. For spatial covariates, run 09_raster_covariates.R and 10_vector_covariates.R\n")
cat("5. Proceed to 01_data_import.R for the site x species analysis\n\n")

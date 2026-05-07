# ============================================================
# 08_moon.R
# Author: Alberto Israel
# Purpose: Compute lunar covariates for each sampling hour
#          using site coordinates and detection timestamps
# Input:   outputs/site_coordinates.csv
#          outputs/bat_data_hourly.csv
# Output:  outputs/lunar_hourly.csv
#          outputs/bat_data_hourly.csv (updated with lunar covariates)
# Notes:
#   - Lunar variables computed via suncalc package
#   - Computed only for hours present in bat_data_hourly.csv
#     (i.e., effectively sampled site × night × hour combinations)
#   - moon_brightness = fraction × moon_visible: ecologically
#     meaningful metric combining illumination and moon above horizon
#   - Aggregation to daily or site level is not performed:
#     lunar conditions vary hourly and averaging would lose
#     biological relevance for bat activity analyses
#   - Run after 00b_data_preparation.R and 07_temperature.R
# Last update: 2026-05-06
# ============================================================


# ------------------------------------------------------------
# 0) Setup project environment
# ------------------------------------------------------------
source(here::here("scripts", "00a_setup.R"))

output_dir <- file.path(here::here(), "outputs")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# ============================================================
# 1) Import input data
# ============================================================

# ------------------------------------------------------------
# 1.1) Site coordinates
# ------------------------------------------------------------

# Produced by 00b_data_preparation.R
# Required columns: SiteID, lat, lon

site_coords <- read_csv(file.path(output_dir, "site_coordinates.csv"),
                        show_col_types = FALSE) |>
  mutate(SiteID = as.integer(SiteID))

cat("Sites loaded:", nrow(site_coords), "\n")
print(site_coords)
cat("\n")

# ------------------------------------------------------------
# 1.2) Hourly detection matrix
# ------------------------------------------------------------

# Produced by 00b_data_preparation.R (updated by 07_temperature.R)
# Required columns: SiteID, Day, Hour

bat_hourly <- read_csv(file.path(output_dir, "bat_data_hourly.csv"),
                       show_col_types = FALSE) |>
  mutate(
    SiteID = as.integer(SiteID),
    Day    = as.Date(Day),
    Hour   = as.integer(Hour)
  )

cat("Hourly detection records loaded:", nrow(bat_hourly), "\n")
cat("Sites:", length(unique(bat_hourly$SiteID)), "\n")
cat("Nights:", length(unique(bat_hourly$Day)), "\n\n")


# ============================================================
# 2) Build sampling grid (SiteID × Day × Hour)
# ============================================================

# Extract unique site × night × hour combinations actually sampled
# This ensures lunar variables are computed only where needed

sampling_grid <- bat_hourly |>
  dplyr::select(SiteID, Day, Hour) |>
  distinct() |>
  left_join(site_coords, by = "SiteID") |>
  arrange(SiteID, Day, Hour)

cat("Unique site × night × hour combinations:", nrow(sampling_grid), "\n\n")

# Check for sites with missing coordinates
missing_coords <- sampling_grid |>
  filter(is.na(lat) | is.na(lon)) |>
  dplyr::select(SiteID) |>
  distinct()

if (nrow(missing_coords) > 0) {
  warning("Missing coordinates for SiteID(s): ",
          paste(missing_coords$SiteID, collapse = ", "),
          "\nThese sites will have NA lunar values.",
          "\nCheck site_coordinates.csv.")
} else {
  cat("All sites have coordinates.\n\n")
}


# ============================================================
# 3) Compute lunar variables
# ============================================================

# suncalc requires POSIXct datetime in UTC
# Construct datetime as the start of each hour (HH:00:00 UTC)

sampling_grid <- sampling_grid |>
  mutate(
    datetime_utc = as.POSIXct(
      paste(Day, sprintf("%02d:00:00", Hour)),
      format = "%Y-%m-%d %H:%M:%S",
      tz = "UTC"
    )
  )

# ------------------------------------------------------------
# 3.1) Moon position (altitude, azimuth, distance)
# ------------------------------------------------------------

# getMoonPosition() requires: date (POSIXct), lat, lon
# Returns altitude and azimuth in radians → convert to degrees

cat("Computing moon position...\n")

moon_position <- getMoonPosition(
  date = sampling_grid$datetime_utc,
  lat  = sampling_grid$lat,
  lon  = sampling_grid$lon
) |>
  as_tibble() |>
  mutate(
    altitude_deg  = round(altitude  * 180 / pi, 4),
    azimuth_deg   = round(azimuth   * 180 / pi, 4),
    distance_km   = round(distance  / 1000,      2)
  ) |>
  dplyr::select(altitude_deg, azimuth_deg, distance_km)

# ------------------------------------------------------------
# 3.2) Moon illumination (fraction, phase, angle)
# ------------------------------------------------------------

# getMoonIllumination() requires: date (POSIXct)
# Returns fraction (0-1), phase (0-1), angle (radians)
# phase: 0 = new moon, 0.25 = first quarter, 0.5 = full, 0.75 = last quarter

cat("Computing moon illumination...\n")

moon_illumination <- getMoonIllumination(
  date = sampling_grid$datetime_utc
) |>
  as_tibble() |>
  mutate(
    fraction = round(fraction, 4),
    phase    = round(phase,    4),
    angle    = round(angle,    4)
  ) |>
  dplyr::select(fraction, phase, angle)

# ------------------------------------------------------------
# 3.3) Combine and derive moon_brightness
# ------------------------------------------------------------

lunar_hourly <- bind_cols(
  sampling_grid |> dplyr::select(SiteID, Day, Hour, datetime_utc),
  moon_position,
  moon_illumination
) |>
  mutate(
    # Moon is visible only when above the horizon (altitude > 0)
    moon_visible    = altitude_deg > 0,
    # moon_brightness: ecologically relevant illumination
    # combines fraction of disk illuminated with visibility above horizon
    # = 0 when moon is below horizon regardless of phase
    moon_brightness = round(fraction * as.integer(moon_visible), 4)
  ) |>
  dplyr::select(SiteID, Day, Hour, datetime_utc,
                altitude_deg, azimuth_deg, distance_km,
                fraction, phase, angle,
                moon_visible, moon_brightness) |>
  arrange(SiteID, Day, Hour)

cat("Lunar variables computed for", nrow(lunar_hourly), "site × night × hour combinations\n\n")


# ============================================================
# 4) Quality checks
# ============================================================

cat("--- Lunar Quality Report ---\n")
cat("Phase range:", round(min(lunar_hourly$phase), 3),
    "to", round(max(lunar_hourly$phase), 3),
    "(0 = new moon, 0.5 = full moon)\n")
cat("Fraction range:", round(min(lunar_hourly$fraction), 3),
    "to", round(max(lunar_hourly$fraction), 3), "\n")
cat("Hours with moon visible:", sum(lunar_hourly$moon_visible),
    "/", nrow(lunar_hourly),
    paste0("(", round(mean(lunar_hourly$moon_visible) * 100, 1), "%)\n"))
cat("Moon brightness range:", round(min(lunar_hourly$moon_brightness), 3),
    "to", round(max(lunar_hourly$moon_brightness), 3), "\n\n")

# Flag any unexpected NA values
n_na <- colSums(is.na(lunar_hourly))
if (any(n_na > 0)) {
  cat("WARNING: NA values detected:\n")
  print(n_na[n_na > 0])
} else {
  cat("No NA values detected.\n\n")
}


# ============================================================
# 5) Save lunar_hourly.csv
# ============================================================

write_csv(lunar_hourly, file.path(output_dir, "lunar_hourly.csv"))
cat("Lunar hourly table saved to outputs/lunar_hourly.csv\n")
cat("Dimensions:", nrow(lunar_hourly), "rows ×", ncol(lunar_hourly), "columns\n\n")


# ============================================================
# 6) Join with hourly detection matrix
# ============================================================

bat_data_hourly <- read_csv(file.path(output_dir, "bat_data_hourly.csv"),
                            show_col_types = FALSE) |>
  mutate(
    SiteID = as.integer(SiteID),
    Day    = as.Date(Day),
    Hour   = as.integer(Hour)
  )

bat_data_hourly <- bat_data_hourly |>
  left_join(
    lunar_hourly |> dplyr::select(SiteID, Day, Hour,
                                   altitude_deg, azimuth_deg, distance_km,
                                   fraction, phase, moon_visible, moon_brightness),
    by = c("SiteID", "Day", "Hour")
  )

# Check join completeness
n_missing <- sum(is.na(bat_data_hourly$moon_brightness))
if (n_missing > 0) {
  cat("WARNING:", n_missing, "rows without lunar match after join\n")
  cat("Check SiteID / Day / Hour alignment between matrices\n\n")
} else {
  cat("Join complete — no missing lunar values.\n\n")
}

write_csv(bat_data_hourly, file.path(output_dir, "bat_data_hourly.csv"))
cat("bat_data_hourly.csv updated with lunar covariates\n\n")


# ============================================================
# 7) Final output summary
# ============================================================

cat("--- Files created/updated in outputs/ ---\n")
cat("  - lunar_hourly.csv          (SiteID × Day × Hour: full lunar variables)\n")
cat("  - bat_data_hourly.csv       (updated with lunar covariates)\n")

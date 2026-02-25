# ============================================================
# 01_data_import.R
# Author: Alberto Israel
# Purpose: Import bat acoustic dataset, clean data, and compute species richness & diversity indices
# Last update: 2026-02-24
# ============================================================


# ------------------------------------------------------------
# 0) Setup project environment
# ------------------------------------------------------------
# Source setup script (installs and loads packages)
source(here::here("scripts", "00_setup.R"))

# Define main directories
data_dir   <- file.path(here::here(), "data")
output_dir <- file.path(here::here(), "outputs")


# ============================================================
# 1) Import data
# ============================================================
bat_data <- read_csv(file.path(data_dir, "bat_data.csv")) #if file is .xlsx, switch to read_excel()
# Inspect
glimpse(bat_data)
summary(bat_data)


# ============================================================
# 2) Clean data
# ============================================================
# Check missing values
bat_data <- bat_data %>%
  filter(!is.na(Site)) %>%    # remove rows with missing site info
  mutate(across(starts_with("species_"), ~replace_na(.x, 0))) # replace NA with 0 for species counts

# Check for duplicates (optional)
bat_data <- bat_data %>%
  distinct()


# ============================================================
# 3) Create community parameters as potential dependet variables
     # for community-level models
# ============================================================

# ------------------------------------------------------------
# Total acoustic activity
# ------------------------------------------------------------
bat_data <- bat_data %>%
  rowwise() %>%
  mutate(
    N_calls = sum(c_across(4:20)), # adjust column range to your dataset (in this case 4 = first species' column, 20 = last species' column)
    N_calls_per_day = N_calls / recording_days) %>% # if recording time is different across samples, standardize to the time unit of interest
  ungroup()


# ------------------------------------------------------------
# Guild-specific acoustic activity (example: open-space foragers)
# ------------------------------------------------------------
bat_data <- bat_data %>%
  mutate(
    N_call_osf = EPTSER + NYCNOC + TADTEN) # Write column namens of OFS; adapt to the OSF species present in your study
# Repeat for edge-space foragers, narrow-space foragers or any other guild designed for your study


# ------------------------------------------------------------
# Species richness
# ------------------------------------------------------------
bat_data <- bat_data %>%
  rowwise() %>%
  mutate(
    N_species = sum(c_across(4:20) > 0) # Adjust column range to your dataset
  ) %>%
  ungroup()


# ------------------------------------------------------------
# Species diversity indices
# ------------------------------------------------------------
bat_data <- bat_data %>%
  mutate(simpson_index = diversity(across(4:20), index = "simpson"), # Adjust column range to your dataset
    shannon_index = diversity(across(4:20), index = "shannon"))

    
# ============================================================
# 4) Quick plots for inspection
# ============================================================


# ------------------------------------------------------------
# Histograms: check overall distribution and potential outliers
# ------------------------------------------------------------
ggplot(bat_data, aes(x = N_call_per_day)) + #switch with species richness or guild-specific number of calls
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  theme_minimal()


# ------------------------------------------------------------
# Boxplots: check distribution and variability of key parameters 
            # per site
# ------------------------------------------------------------

variables_to_plot <- c("N_call_per_day", "N_species", "shannon_index", "simpson_index")

for (var in variables_to_plot) {
  ggplot(bat_data, aes_string(x = "Site", y = var)) +
    geom_boxplot(fill = "lightgreen") +
    theme_minimal() +
    labs(title = paste("Distribution of", var, "per site"), x = "Site", y = var)
}


# ------------------------------------------------------------
# Scatterplots: check relationships between key variables
# ------------------------------------------------------------

# Example: Guild-specific activity vs Total activity
ggplot(bat_data, aes(x = N_calls, y = N_call_osf)) +
  geom_point(color = "darkorange") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal()


# ============================================================
# 5) Save cleaned dataset for analysis
# ============================================================
write_csv(bat_data, file.path(output_dir, "bat_data_cleaned.csv"))

# ============================================================
# 05_additional_analyses.R
# Author: Alberto Israel
# Purpose: Additional fine-scale community analyses
# Section a: Species co-occurrence analysis (habitat-stratified)
# Notes:
#   - Requires a site × species matrix
#   - Presence/absence-based co-occurrence
#   - Designed to complement 03 (univariate models) and 04 (community-level inference)
# Last update: 2026-02-24
# ============================================================

# ------------------------------------------------------------
# 0) Setup project environment
# ------------------------------------------------------------
source(here::here("scripts", "00_setup.R"))

data_dir   <- file.path(here::here(), "data")
output_dir <- file.path(here::here(), "outputs")

library(tidyverse)  # data manipulation and plotting
library(cooccur)    # co-occurrence analysis
library(grid)       # text annotation on plots


# ------------------------------------------------------------
# 1) Load prepared dataset
# ------------------------------------------------------------
bat_data <- read_csv(file.path(output_dir, "bat_data_cleaned.csv"))

# Column range of species matrix (site × species)
species_cols <- 4:20   # <<--- FILL: first:last species column

# Site identifier
site_id <- "site"      # <<--- FILL

# Habitat / grouping variable for stratification
habitat_var <- "habitat_class"   # <<--- FILL (e.g. "Forest" vs "Open")


# Pre-analysis column check
required_cols <- c(site_id, habitat_var)
missing_cols <- setdiff(required_cols, colnames(bat_data))
if(length(missing_cols) > 0) stop(paste("Missing columns in dataset:", paste(missing_cols, collapse=", ")))


# ------------------------------------------------------------
# 2) Helper function: build presence/absence matrix for cooccur
# ------------------------------------------------------------

prepare_cooccur_matrix <- function(df, species_cols, site_id) {
  
  mat <- as.matrix(df[, species_cols])
  rownames(mat) <- df[[site_id]]
  
  # Convert abundance to presence/absence
  mat <- ifelse(mat > 0, 1, 0)
  
  # Transpose: species (rows) × sites (columns)
  mat <- t(mat)
  
  # Remove species with zero occurrences
  mat <- mat[rowSums(mat) > 0, ]
  
  # Check if matrix is empty
  if(nrow(mat) == 0 | ncol(mat) == 0) stop("Co-occurrence matrix empty after filtering. Check species_cols and site_id")
  
  return(mat)
}

# ------------------------------------------------------------
# 3) Co-occurrence analysis by habitat (stratified)
# ------------------------------------------------------------

run_cooccurrence_by_group <- function(df, group_value, habitat_var, species_cols, site_id) {
  
  df_sub <- df %>% 
    filter(.data[[habitat_var]] == group_value)
  
  mat_sub <- prepare_cooccur_matrix(df_sub, species_cols, site_id)
  
  coc_result <- cooccur(mat_sub, spp_names = TRUE, thresh = FALSE)
  
  return(coc_result)
}

# Identify habitat groups
habitat_levels <- na.omit(unique(bat_data[[habitat_var]]))

# Run co-occurrence for each habitat
cooccur_results <- lapply(habitat_levels, function(h) {
  df_sub <- bat_data %>% filter(.data[[habitat_var]] == h)
  
  if(nrow(df_sub) == 0) {
    warning(paste("Habitat", h, "empty, skipped"))
    return(NULL)
  }
  
  run_cooccurrence_by_group(df_sub, h, habitat_var, species_cols, site_id)
})

names(cooccur_results) <- habitat_levels

# ------------------------------------------------------------
# 4) Summaries and plots (per habitat)
# ------------------------------------------------------------

for (h in habitat_levels) {
  
  coc <- cooccur_results[[h]]
  
  if(!is.null(coc)) {
    cat("\n========================================\n")
    cat("Habitat:", h, "\n")
    cat("========================================\n")
    
    print(summary(coc))
    
    plot(coc)
    grid.text(h, x = 0.05, y = 0.95, just = c("left", "top"),
              gp = gpar(fontsize = 14, fontface = "bold"))
  } else {
    warning(paste("Habitat", h, "is NULL and was skipped in plotting"))
  }
}

# ------------------------------------------------------------
# 5) Co-occurrence on full dataset (reference)
# ------------------------------------------------------------

mat_all <- prepare_cooccur_matrix(bat_data, species_cols, site_id)

if(ncol(mat_all) == 0) stop("Empty matrix: check species_cols and site_id")

if(nrow(mat_all) == 0) stop("Empty matrix: no sites after filtering")

coc_all <- cooccur(mat_all, spp_names = TRUE, thresh = FALSE)

cat("\n========================================\n")
cat("All habitats combined\n")
cat("========================================\n")

print(summary(coc_all))
plot(coc_all)
grid.text("All habitats", x = 0.05, y = 0.95, just = c("left", "top"),
          gp = gpar(fontsize = 14, fontface = "bold"))


# ------------------------------------------------------------
# 6) Extract results tables (optional: for comparison or export)
# ------------------------------------------------------------

for(h in habitat_levels) {
  if(!is.null(cooccur_results[[h]])) {
    cooccur_tables <- cooccur_results[[h]]$results
    write_csv(cooccur_tables, file.path(output_dir, paste0("cooccurrence_", h, ".csv")))
    message(paste("Saved co-occurrence table for habitat:", h))
  }
}


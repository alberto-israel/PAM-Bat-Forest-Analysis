# ============================================================
# 05_additional_analyses.R
# Author: Alberto Israel
# Purpose: Additional fine-scale community analyses
# Section b: : Beta diversity (turnover, nestedness, total, PERMANOVA)
# Notes:
#   - Designed to integrate with 01_data_import.R
# Last update: 2026-02-24
# ============================================================

# -----------------------------
# 0) Setup project environment
# -----------------------------
source(here::here("scripts", "00_setup.R"))

library(tidyverse)  # data manipulation and plotting
library(vegan)      # PERMANOVA, ordination
library(betapart)   # beta diversity partitioning

data_dir   <- file.path(here::here(), "data")
output_dir <- file.path(here::here(), "outputs")

# -----------------------------
# 1) # Load prepared dataset
# -----------------------------
bat_data <- read_csv(file.path(output_dir, "bat_data_cleaned.csv"))

# Column range of species matrix (site × species)
species_cols <- 4:19  # <<--- FILL: first:last species column

# Environmental / predictor variables
env_vars <- c("predictor1", "predictor2")  # <<--- FILL: numeric or factor predictors

# Optional random factors for mixed models
random_factors <- c("site", "night")  # <<--- FILL if needed

# -----------------------------
# 2) Prepare data
# -----------------------------
# Check columns
required_species_cols <- species_cols
if(any(required_species_cols > ncol(bat_data))) stop("species_cols exceed dataset columns")

if(!all(env_vars %in% colnames(bat_data))) stop("Some env_vars not found in dataset")

# Species matrix (site × species)
species_matrix <- bat_data[, species_cols]

# Convert to presence/absence (0/1), handle NAs
species_matrix[species_matrix > 0] <- 1
species_matrix[is.na(species_matrix)] <- 0

# Row names: site IDs
site_id <- "site"  # <<--- FILL
if(!site_id %in% colnames(bat_data)) stop("site_id column not found in dataset")
rownames(species_matrix) <- bat_data[[site_id]]

# Environmental matrix
env_data <- bat_data[, env_vars]

# -----------------------------
# 3) Compute beta diversity
# -----------------------------
if(nrow(species_matrix) == 0 | ncol(species_matrix) == 0) stop("Empty species matrix")

# Partition beta diversity into turnover and nestedness
beta_res <- beta.pair(species_matrix, index.family = "sorensen") # Sorensen-based
turnover <- beta_res$beta.sim
nestedness <- beta_res$beta.sne
total_beta <- beta_res$beta.sor

# -----------------------------
# 4) Summarize
# -----------------------------
cat("Mean turnover:", mean(as.vector(turnover), na.rm = TRUE), "\n")
cat("Mean nestedness:", mean(as.vector(nestedness), na.rm = TRUE), "\n")
cat("Mean total beta:", mean(as.vector(total_beta), na.rm = TRUE), "\n")

# -----------------------------
# 5) Optional: visualize dissimilarity matrices
# -----------------------------
# Example heatmap for total beta
heatmap(as.matrix(total_beta), symm = TRUE, main = "Total Beta Diversity (Sorensen)") 

# -----------------------------
# 6) Relate beta diversity to environmental variables
# -----------------------------
# PERMANOVA on total beta
adonis_res <- adonis2(as.dist(total_beta) ~ ., data = env_data, permutations = 999)
print(adonis_res)

# -----------------------------
# 7) Notes
# -----------------------------
# - species_matrix: site x species, presence/absence (0/1)
# - env_data: numeric or factor variables
# - For turnover/nestedness visualizations, consider using heatmaps or ordinations
# - Replace column indices and names with your dataset-specific columns

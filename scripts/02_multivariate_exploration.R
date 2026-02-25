# ============================================================
# 02_multivariate_exploration.R
# Author: Alberto Israel
# Purpose: Exploratory multivariate analysis of bat community composition
# Methods: PCA (Hellinger-transformed data) and NMDS (Bray–Curtis)
# Notes:
#   - Unconstrained ordination only (no formal modelling)
#   - Designed to integrate with 01_data_import.R
# Last update: 2026-02-24
# ============================================================


# ------------------------------------------------------------
# 0) Setup project environment
# ------------------------------------------------------------
source(here::here("scripts", "00_setup.R"))

data_dir   <- file.path(here::here(), "data")
output_dir <- file.path(here::here(), "outputs")


# ============================================================
# 1) Load prepared dataset
# ============================================================
bat_data <- read_csv(file.path(output_dir, "bat_data_cleaned.csv"))

# Check
glimpse(bat_data)


# ============================================================
# 2) Define community (species) matrix
# ============================================================
# NOTE:
# Species data are assumed to be stored in contiguous columns.
# Modify 'species_cols' to match the first and last species columns
# in your dataset (e.g. 4:20).

species_cols <- 4:20
species_matrix <- bat_data[, species_cols]

# Check structure
summary(species_matrix)

# Optional: remove rows with all-zero species counts
species_matrix <- species_matrix[rowSums(species_matrix, na.rm = TRUE) > 0, ]


# ============================================================
# 3) Data transformation for ordination
# ============================================================
# Hellinger transformation is recommended for Euclidean-based ordinations (e.g. PCA)
species_hellinger <- decostand(species_matrix, method = "hellinger")


# ============================================================
# 4) PCA: Exploratory ordination (Hellinger-transformed data)
# ============================================================

pca_result <- dudi.pca(species_hellinger, scannf = FALSE, nf = 2)

# Basic scatterplot
scatter(pca_result, clab.row = 0,
        main = "PCA of Bat Community Composition (Hellinger-transformed data)")

# Biplot with factoextra
fviz_pca_biplot(
  pca_result,
  repel = TRUE,
  title = "PCA Biplot")


# ------------------------------------------------------------
# 4a) Example grouping by a categorical environmental variable
# ------------------------------------------------------------
# NOTE:
# This block is optional and purely exploratory.
# Replace 'env_class' with any categorical environmental variable available
# in your dataset (e.g. habitat type, management class, forest type).


if ("env_class" %in% colnames(bat_data)) {
  
  s.class(
    pca_result$li,
    fac = factor(bat_data$env_class),
    col = rainbow(length(unique(bat_data$env_class))),
    sub = "Sites grouped by environmental category (example variable)"
  )
} else {
  message("No categorical environmental variable found (e.g. 'env_class'). Skipping group plot.")
}


# ============================================================
# 5) NMDS: Community dissimilarity (Bray–Curtis)
# ============================================================

# Compute Bray–Curtis dissimilarity on raw abundance data
bray_dist <- vegdist(species_matrix, method = "bray")

# Run NMDS (2 dimensions) with sufficient iterations and reproducibility
set.seed(123)
nmds_result <- metaMDS(bray_dist, k = 2, trymax = 100, autotransform = FALSE)

# Check stress value
nmds_result$stress

# Plot NMDS
plot(nmds_result, type = "t",
     main = "NMDS of Bat Community Composition (Bray–Curtis)")


# ------------------------------------------------------------
# 5a) Fit environmental variables onto NMDS (optional, exploratory)
# ------------------------------------------------------------
# NOTE:
# This block is optional and purely exploratory.
# Replace 'env_class' with any categorical environmental variable available
# in your dataset (e.g. habitat type, management class, forest type)

if ("env_class" %in% colnames(bat_data)) {
  
  env_data <- bat_data %>%
    select(env_class)
  
  envfit_result <- envfit(nmds_result, env_data, permutations = 999)
  
  plot(envfit_result, p.max = 0.05, col = "black")
} else {
  message("No environmental variables available for envfit(). Skipping.")
}


# ------------------------------------------------------------
# 5b) Mean dissimilarity (descriptive)
# ------------------------------------------------------------
mean_dissimilarity <- mean(bray_dist)
mean_dissimilarity


# ============================================================
# 6) Additional descriptive summaries of community structure
# ============================================================

# Total number of calls per species
total_calls_per_species <- colSums(species_matrix, na.rm = TRUE)

# Number of sites where each species is present
site_occurrence_per_species <- sapply(
  species_matrix,
  function(x) sum(x > 0, na.rm = TRUE)
)

# Inspect summaries
total_calls_per_species
site_occurrence_per_species


# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

species_summary <- data.frame(
  species = names(total_calls_per_species),
  total_calls = total_calls_per_species,
  site_occurrence = site_occurrence_per_species
)

write_csv(species_summary, file.path(output_dir, "species_summary.csv"))

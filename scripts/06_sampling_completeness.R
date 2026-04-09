# ============================================================
# 06_sampling_completeness.R
# Author: Alberto Israel
# Purpose: Evaluate sampling completeness and adequacy for bat PAM data
# Methods: Species accumulation curves (spatial and temporal),
#          richness estimators (Chao1, Jackknife), coverage-based 
#          rarefaction (iNEXT), stratified completeness analysis
# Notes:
#   - Temporal accumulation (Section 4) requires daily-level data
#   - If daily data unavailable, Section 4 will be skipped
#   - Designed to integrate with 01_data_import.R
# Last update: 2026-03-20
# ============================================================

# ------------------------------------------------------------
# 0) Setup project environment
# ------------------------------------------------------------
source(here::here("scripts", "00_setup.R"))

# Additional package for coverage-based rarefaction
if (!require("iNEXT")) install.packages("iNEXT")
library(iNEXT)

data_dir   <- file.path(here::here(), "data")
output_dir <- file.path(here::here(), "outputs")

# ------------------------------------------------------------
# 0.1) Load prepared dataset
# ------------------------------------------------------------

bat_data <- read_csv(file.path(output_dir, "bat_data_cleaned.csv"))

# Check
glimpse(bat_data)

# ------------------------------------------------------------
# 0.2) Define community (species) matrix
# ------------------------------------------------------------

# NOTE:
# Species data are assumed to be stored in contiguous columns.
# Modify 'species_cols' to match the first and last species columns
# in your dataset (e.g. 4:20).

species_cols <- 4:20.    # adapt to your dataset
species_matrix <- bat_data[, species_cols]

# Check structure
summary(species_matrix)

# Optional: remove rows with all-zero species counts
species_matrix <- species_matrix[rowSums(species_matrix, na.rm = TRUE) > 0, ]

# Create presence/absence matrix for accumulation analyses
species_pa <- ifelse(species_matrix > 0, 1, 0)

# ============================================================
#1) Species Accumulation Curve - Sites
# ============================================================

# Compute species accumulation curve with random method
spec_accum_sites <- specaccum(species_pa, method = "random", permutations = 999)

# Plot
plot(spec_accum_sites, 
     xlab = "Number of sampling sites", 
     ylab = "Cumulative species richness",
     ci.type = "poly", col = "black", lwd = 2,
     ci.lty = 0, ci.col = "lightgrey",
     main = "Species Accumulation Curve - Spatial Sampling")

# Check if plateau is reached
final_richness <- tail(spec_accum_sites$richness, 1)
final_se <- tail(spec_accum_sites$sd, 1)

# ============================================================
#2) Richness Estimators
# ============================================================

# Compute non-parametric richness estimators
estimators <- specpool(species_pa)
print(estimators)

# Extract key values
obs_rich <- sum(colSums(species_matrix) > 0)
chao_estimate <- estimators$chao
chao_se <- estimators$chao.se
jack1_estimate <- estimators$jack1
jack1_se <- estimators$jack1.se

# Calculate sampling completeness
completeness_chao <- (obs_rich / chao_estimate) * 100
completeness_jack1 <- (obs_rich / jack1_estimate) * 100

# Summary
cat("\n--- Sampling Completeness Summary ---\n")
cat("Observed richness:       ", obs_rich, "\n")
cat("Chao1 estimate:          ", round(chao_estimate, 1), " ± ", round(chao_se, 1), " SE\n", sep = "")
cat("Jack1 estimate:          ", round(jack1_estimate, 1), " ± ", round(jack1_se, 1), " SE\n", sep = "")
cat("Completeness (Chao1):    ", round(completeness_chao, 1), "%\n", sep = "")
cat("Completeness (Jack1):    ", round(completeness_jack1, 1), "%\n\n", sep = "")

# Visualization: Observed vs Estimated
barplot(c(Observed = obs_rich, 
          Chao1 = chao_estimate,
          Jack1 = jack1_estimate),
        col = c("black", "grey60", "grey80"),
        ylab = "Species richness",
        main = "Observed vs Estimated Species Richness",
        ylim = c(0, max(jack1_estimate) * 1.1))
abline(h = obs_rich, lty = 2, col = "red", lwd = 1.5)
legend("topright", 
       legend = paste0("Completeness (Chao1): ", round(completeness_chao, 1), "%"),
       bty = "n", cex = 0.9)

# ============================================================
#3) Coverage-based Rarefaction (iNEXT)
# ============================================================

# Prepare data for iNEXT (transpose to species × sites)
data_inext <- t(species_pa)

# Run iNEXT for species richness (q = 0)
inext_result <- iNEXT(data_inext, 
                      q = 0,
                      datatype = "incidence_raw",
                      endpoint = 2 * ncol(data_inext))

# Plot rarefaction and extrapolation curves with confidence intervals
ggiNEXT(inext_result, type = 1) +
  theme_bw() +
  labs(title = "Sample Coverage and Species Richness",
       x = "Number of sampling units",
       y = "Species richness") +
  theme(legend.position = "bottom")

# Extract and report sample coverage
coverage_info <- inext_result$iNextEst$coverage_based
sample_coverage <- coverage_info$SC[coverage_info$Order.q == 0 & 
                                      coverage_info$Method == "Observed"][1]

cat("Sample coverage (observed): ", round(sample_coverage, 3), "\n\n", sep = "")

# ============================================================
#4) Temporal Accumulation - Days (OPTIONAL)
# ============================================================
# NOTE:
# This analysis requires day-level data (site × day × species).
# If you only have aggregated site × species data, this section will be skipped.

# Check for daily data file
daily_file <- file.path(data_dir, "bat_data_daily.csv")

if (file.exists(daily_file)) {     
  
  message("Daily data found. Running temporal accumulation analysis...\n")
  
# ------------------------------------------------------------
# 4a) Load and prepare daily data
# ------------------------------------------------------------
  
  data_perday <- read_csv(daily_file)
  
  # NOTE: Adjust column indices to match your daily dataset structure
  # Expected format: Site | Day | Species1 | Species2 | ...
  spec_cols_daily <- 3:ncol(data_perday)  # Adjust if needed
  
  # Convert to presence/absence
  dati_pa_perday <- data_perday %>%
    mutate(across(all_of(spec_cols_daily), ~ifelse(. > 0, 1, 0)))
  
# ------------------------------------------------------------
# 4b) Compute accumulation curves for each site
# ------------------------------------------------------------
  
  siti <- unique(dati_pa_perday$Site)
  
  curves_per_site <- lapply(siti, function(s) {
    tmp <- dati_pa_perday %>% 
      filter(Site == s) %>% 
      arrange(Day)
    
    mat <- as.matrix(tmp[, spec_cols_daily])
    specaccum(mat, method = "collector")
  })
  
# ------------------------------------------------------------
# 4c) Aggregate across sites: mean richness by day
# ------------------------------------------------------------
  
  # Find maximum number of sampling days across all sites
  max_days <- max(sapply(curves_per_site, function(x) length(x$richness)))
  
  # Build matrix: rows = sites, columns = days
  richness_matrix <- sapply(curves_per_site, function(x){
    r <- x$richness
    # If curve is shorter than max_days, pad with last value
    if(length(r) < max_days) {
      r <- c(r, rep(r[length(r)], max_days - length(r)))
    }
    return(r)
  })
  
  richness_matrix <- t(richness_matrix)
  
  # Calculate mean and standard error across sites
  mean_richness <- apply(richness_matrix, 2, mean)
  sd_richness <- apply(richness_matrix, 2, sd)
  se_richness <- sd_richness / sqrt(nrow(richness_matrix))
  
  days <- 1:max_days
  
  # ------------------------------------------------------------
  # 4d) Plot temporal accumulation curve
  # ------------------------------------------------------------
  
  # Plot with custom axes (integer ticks only)
  plot(days, mean_richness, type = "l", lwd = 2, col = "black",
       xlab = "Recording days per site",
       ylab = "Cumulative species richness",
       main = "Mean Species Accumulation Curve - Temporal Sampling",
       xaxt = "n", yaxt = "n")
  
  # Add ±SE envelope
  lines(days, mean_richness + se_richness, lty = 2, col = "black")
  lines(days, mean_richness - se_richness, lty = 2, col = "black")
  
  # Custom integer axes
  axis(1, at = days, labels = days)
  y_ticks <- floor(min(mean_richness - se_richness)):ceiling(max(mean_richness + se_richness))
  axis(2, at = y_ticks, labels = y_ticks)
  
  # Report summary
  cat("Number of sites analyzed: ", length(siti), "\n")
  cat("Mean richness after 1 day:  ", round(mean_richness[1], 1), " ± ", round(se_richness[1], 2), " SE\n", sep = "")
  cat("Mean richness at max days (", max_days, "): ", round(mean_richness[max_days], 1), " ± ", round(se_richness[max_days], 2), " SE\n\n", sep = "")
  
} else {
  
  message("⚠️  Daily data not found: ", daily_file)
  message("Skipping temporal accumulation analysis.")
  message("\nTo run this analysis:")
  message("  1. Prepare a dataset with columns: Site, Day, Species1, Species2, ...")
  message("  2. Save as 'bat_data_daily.csv' in the data/ folder")
  message("  3. Re-run this script\n")
  
}

# ============================================================
# SECTION 5: Stratified Completeness - by Habitat
# ============================================================
# NOTE:
# This section assumes a 'habitat' or similar categorical variable exists.
# Adjust the variable name to match your dataset (e.g. 'forest_type', 'land_cover').

if ("habitat" %in% colnames(bat_data)) {
  
  habitats <- unique(bat_data$habitat)
  completeness_habitat <- data.frame()
  
  for (hab in habitats) {
    
    # Subset data for current habitat
    subset_idx <- bat_data$habitat == hab
    subset_pa <- species_pa[subset_idx, ]
    
    # Compute estimators
    est <- specpool(subset_pa)
    obs <- sum(colSums(subset_pa) > 0)
    comp_chao <- (obs / est$chao) * 100
    
    # Store results
    completeness_habitat <- rbind(completeness_habitat,
                                  data.frame(
                                    Habitat = hab,
                                    N_sites = sum(subset_idx),
                                    Observed = obs,
                                    Chao1 = round(est$chao, 1),
                                    Completeness = round(comp_chao, 1)
                                  ))
  }
  
  print(completeness_habitat)
  
  # Visualization
  ggplot(completeness_habitat, aes(x = Habitat, y = Completeness)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = 90, linetype = "dashed", color = "red") +
    geom_text(aes(label = paste0(round(Completeness, 1), "%")), 
              vjust = -0.5, size = 4) +
    ylim(0, 105) +
    theme_minimal() +
    labs(y = "Sampling completeness (%)",
         x = "Habitat type",
         title = "Sampling Completeness by Habitat",
         caption = "Red dashed line = 90% completeness threshold") +
    theme(plot.caption = element_text(hjust = 0))
  
} else {
  
  message("No 'habitat' variable found in the dataset.")
  message("Skipping stratified analysis.")
  message("To run this analysis, ensure your dataset includes a categorical habitat variable.\n")
  
}

# ============================================================
#6) Summary and Recommendations
# ============================================================


# Compile summary table
summary_completeness <- data.frame(
  Metric = c("Observed richness", 
             "Chao1 estimate", 
             "Jack1 estimate",
             "Completeness (Chao1, %)",
             "Completeness (Jack1, %)",
             "Sample coverage"),
  Value = c(obs_rich, 
            round(chao_estimate, 1),
            round(jack1_estimate, 1),
            round(completeness_chao, 1),
            round(completeness_jack1, 1),
            round(sample_coverage, 3))
)

print(summary_completeness)

# Save summary
write_csv(summary_completeness, 
          file.path(output_dir, "sampling_completeness_summary.csv"))

message("\nSummary saved to: ", file.path(output_dir, "sampling_completeness_summary.csv"))

# ------------------------------------------------------------
# Interpretation and recommendations
# ------------------------------------------------------------

cat("\n--- Interpretation ---\n")

if (completeness_chao < 80) {
  cat("Warning: Sampling completeness is below 80% (", round(completeness_chao, 1), "%).\n", sep = "")
  cat("   Recommendation: Consider additional sampling effort to improve species inventory.\n")
  cat("   Estimated missing species: ", round(chao_estimate - obs_rich, 1), "\n\n", sep = "")
  
} else if (completeness_chao < 90) {
  cat("NOTE: Sampling completeness is adequate (", round(completeness_chao, 1), "%) but could be improved.\n", sep = "")
  cat("   Recommendation: Additional sampling may detect ", round(chao_estimate - obs_rich, 1), " more species.\n\n", sep = "")
  
} else {
  cat("Positive: Sampling completeness exceeds 90% (", round(completeness_chao, 1), "%).\n", sep = "")
  cat("   Sampling effort is adequate for characterizing the bat community.\n\n")
}

cat("Sample coverage: ", round(sample_coverage, 3), "\n", sep = "")
if (sample_coverage < 0.95) {
  cat(" Low coverage suggests rare species are likely under-represented.\n\n")
} else {
  cat("High coverage indicates sampling captured most species in the assemblage.\n\n")
}

message("=== Analysis complete ===\n")

# ============================================================
# 04_multivariate_inference.R
# Author: Alberto Israel
# Purpose: Inferential multivariate analysis of bat community composition
# Methods: PERMANOVA, RDA / CCA, Variance Partitioning, Indicator Species Analysis
# Notes:
#   - Designed to integrate with:
#       01_data_import.R
#       02_multivariate_exploration.R
#       03_modeling_univariate.R
#   - Fill in the placeholders for response and predictors
# Last update: 2026-02-24
# ============================================================


# ------------------------------------------------------------
# 0) Setup project environment
# ------------------------------------------------------------
source(here::here("scripts", "00_setup.R"))

data_dir   <- file.path(here::here(), "data")
output_dir <- file.path(here::here(), "outputs")


# ------------------------------------------------------------
# 0) Load prepared dataset
# ------------------------------------------------------------
bat_data <- read_csv(file.path(output_dir, "bat_data_cleaned.csv"))


# ------------------------------------------------------------
# 0) Define community matrix (SPECIFY COLUMN RANGE)
# ------------------------------------------------------------
species_columns <- 4:20 # Replace column indices with the correct range for species data
species_matrix <- bat_data[, species_columns]

# Check
summary(species_matrix)


# ------------------------------------------------------------
# 0) Define environmental predictors (FILL THESE)
# ------------------------------------------------------------
# Example placeholders – replace with your actual predictors
predictors <- c("your", "list", "of_predictors")   # e.g. c("canopy_cover", "edge_density", "ndvi")

env_data <- bat_data %>%
  select(all_of(predictors))


# ------------------------------------------------------------
# 0) Data transformation & assumptions
# ------------------------------------------------------------

# Hellinger transformation (recommended for community data)
species_hellinger <- decostand(species_matrix, method = "hellinger")

# Check multicollinearity among predictors
if (ncol(env_data) > 1) {
  temp_formula <- as.formula(paste("~", paste(names(env_data), collapse = " + ")))
  vif_values <- vif(lm(temp_formula, data = env_data))
  print(vif_values)
} else {
  message("VIF requires at least 2 predictors. Skipping.")
}


# ------------------------------------------------------------
# 0) Gradient length: DCA to choose RDA vs CCA
# ------------------------------------------------------------
dca_result <- decorana(species_hellinger)
dca_result

# Rule of thumb:
#   - Axis length < 3 → RDA (linear)
#   - Axis length > 4 → CCA (unimodal)



# ============================================================
# 1) PERMANOVA: Community-level hypothesis testing
# ============================================================

# 1.1 Distance matrix
bray_dist <- vegdist(species_matrix, method = "bray")

# 1.2 PERMANOVA model
permanova_result <- adonis2(bray_dist ~ ., data = env_data, permutations = 999)
permanova_result

# 1.3 Check homogeneity of dispersions
dispersion_test <- betadisper(bray_dist, groups = env_data[[predictors[1]]])
anova(dispersion_test)


# ============================================================
# 2) Constrained ordination: RDA or CCA
# ============================================================

# 2.1 Choose method based on DCA results
# If gradient length < 3 → use RDA
# If gradient length > 4 → use CCA
grad_len <- max(dca_result$evals)
if (grad_len < 3) {
  ordination_method <- "RDA"
  message("DCA axis length < 3: using RDA (linear gradients)")
} else if (grad_len > 4) {
  ordination_method <- "CCA"
  message("DCA axis length > 4: using CCA (unimodal gradients)")
} else {
  ordination_method <- "RDA"  # default for 3-4 zone
  message("DCA axis length 3-4 (intermediate): defaulting to RDA")
}


if (ordination_method == "RDA") {
  
  ord_model <- rda(species_hellinger ~ ., data = env_data)
  
} else if (ordination_method == "CCA") {
  
  ord_model <- cca(species_matrix ~ ., data = env_data)
  
}

# 2.2 Global test
anova_global <- anova(ord_model, permutations = 999)
anova_global

# 2.3 Test by axis
anova_axis <- anova(ord_model, by = "axis", permutations = 999)
anova_axis

# 2.4 Test by terms (predictors)
anova_terms <- anova(ord_model, by = "terms", permutations = 999)
anova_terms

# 2.5 Ordination plot
plot(ord_model, scaling = 2, main = paste(ordination_method, "of Bat Community"))



# ============================================================
# 3) Variance partitioning (optional: if multiple predictor sets)
# ============================================================

# Example: define predictor blocks (EDIT OR REMOVE IF NOT NEEDED)
env_local <- bat_data %>% select(your_local_variables)
env_landscape <- bat_data %>% select(your_landscape_variables)

varpart_result <- varpart(
  species_hellinger,
  env_local,
  env_landscape)

varpart_result
plot(varpart_result, bg = c("lightblue", "lightgreen"), Xnames = c("Local", "Landscape"))



# ============================================================
# 4) Indicator Species Analysis
# ============================================================

# 4.1 Define grouping variable (FILL)
grouping_variable <- "env_class"   # e.g. habitat type, management class, forest type

if (grouping_variable %in% colnames(bat_data)) {
  
  groups <- factor(bat_data[[grouping_variable]])
  
  # 4.2 Run IndVal
  indval_result <- multipatt(
    species_matrix,
    groups,
    func = "IndVal.g",
    control = how(nperm = 999))
  
  # 4.3 Summary of significant indicator species
  summary(indval_result, alpha = 0.05, indvalcomp = TRUE)
  
} else {
  message("Grouping variable not found. Indicator species analysis skipped.")
}



# ============================================================
# 5) Sensitivity checks (optional but recommended)
# ============================================================

# 5.1 Presence / absence transformation
species_pa <- decostand(species_matrix, method = "pa")

bray_pa <- vegdist(species_pa, method = "bray")

permanova_pa <- adonis2(
  bray_pa ~ .,
  data = env_data,
  permutations = 999)
permanova_pa


# ============================================================
# 6) Save key outputs
# ============================================================

# 6.1 Save PERMANOVA table
write_csv(as.data.frame(permanova_result), file.path(output_dir, "permanova_results.csv"))

# 6.2 Save ordination site scores
site_scores <- scores(ord_model, display = "sites", scaling = 2) %>% as.data.frame()
write_csv(site_scores, file.path(output_dir, "ordination_site_scores.csv"))

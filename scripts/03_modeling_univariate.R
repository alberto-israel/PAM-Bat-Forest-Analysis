# ============================================================
# 03_modeling_univariate.R
# Author: Alberto Israel
# Purpose: Univariate and species-level modeling of community metrics
# Methods: GLMs, GLMMs and GAMs
# Notes:
#   - Designed to integrate with 01_data_import.R and 02_multivariate_exploration.R
#   - All variable names must be adapted to the user's dataset
# Last update: 2026-02-24
# ============================================================

# ------------------------------------------------------------
# 0.0) Setup project environment
# ------------------------------------------------------------
source(here::here("scripts", "00_setup.R"))

data_dir   <- file.path(here::here(), "data")
output_dir <- file.path(here::here(), "outputs")

# ============================================================
# 0.1) Load prepared dataset
# ============================================================
bat_data <- read_csv(file.path(output_dir, "bat_data_cleaned.csv"))


# ----------------------------------------------------------
# USER-DEFINED VARIABLES (TO FILL)
# ----------------------------------------------------------

response_activity <- "Your_activity_variable"    # e.g. total calls
response_richness <- "Your_richness_variable"    # e.g. species richness
response_evenness <- "Your_evenness_variable"    # e.g. Simpson or Shannon

predictors <- c("your", "list", "of_predictors") # e.g. habitat, structure, landscape

random_site  <- "site"   # optional random effect for GLMMs
random_night <- "night"  # optional random effect for GLMMs

# Optional: species columns for species-level analysis
species_columns <- c("SPECIES_1", "SPECIES_2", "SPECIES_3")  # fill with your species columns


# ==========================================================
# 01) COLLINEARITY AND EXPLORATORY CHECKS
# ==========================================================

# 1.1) Correlation matrix (Spearman, robust to non-normality)
cor_mat <- cor(bat_data[, predictors], method = "spearman", use = "pairwise.complete.obs")
corrplot(cor_mat, method = "color", tl.col = "black", addCoef.col = "black")

# 1.2) Variance Inflation Factor (true multicollinearity)
vif_lm <- lm(as.formula(paste(response_activity, "~", paste(predictors, collapse = " + "))),
             data = bat_data)
check_collinearity(vif_lm)


# ==========================================================
# 02) TOTAL ACTIVITY MODELING
# ==========================================================


# ------------------------------------------------------------
# 2.1) Poisson GLM
# ------------------------------------------------------------
m_act_pois <- glm(
  as.formula(paste(response_activity, "~", paste(predictors, collapse = " + "))),
  family = poisson, data = bat_data)
summary(m_act_pois)
check_overdispersion(m_act_pois)


# ------------------------------------------------------------
# 2.2) Negative Binomial GLM
# ------------------------------------------------------------
m_act_nb <- glm.nb(
  as.formula(paste(response_activity, "~", paste(predictors, collapse = " + "))),
  data = bat_data)
summary(m_act_nb)

check_overdispersion(m_act_nb)


# ------------------------------------------------------------
# 2.3) GLMM (if random effects are available)
# ------------------------------------------------------------
if (random_site %in% colnames(bat_data)) {
  m_act_glmm <- glmmTMB(
    as.formula(paste(response_activity, "~",
                     paste(predictors, collapse = " + "),
                     "+ (1 |", random_site, ")")),
    family = nbinom2, data = bat_data)
  summary(m_act_glmm)
  AIC(m_act_pois, m_act_nb, m_act_glmm)
} else {
  AIC(m_act_pois, m_act_nb)
}


# ------------------------------------------------------------
# 2.4) Diagnostics
# ------------------------------------------------------------
par(mfrow = c(2,2))
plot(m_act_nb)
which(cooks.distance(m_act_nb) > 4 / nrow(bat_data))

res.model <- simulateResiduals(m_act_glmm)
plot(res.model)
testDispersion(res.model)
testZeroInflation(res.model)
testResiduals(res.model)


# ------------------------------------------------------------
# 2.5) Visualization
# ------------------------------------------------------------
plot_model(m_act_nb, type = "est", transform = NULL, show.values = TRUE, show.p = TRUE) +
  theme_bw() + labs(title = "Total activity model", y = "Predictors", x = "Estimates")

lapply(predictors, function(v)
  visreg(m_act_nb, v, scale = "response",
         xlab = v, ylab = response_activity))



# ==========================================================
# 03) SPECIES RICHNESS MODELING
# ==========================================================


# ------------------------------------------------------------
# 3.1) Poisson GLM
# ------------------------------------------------------------
m_rich_pois <- glm(
  as.formula(paste(response_richness, "~", paste(predictors, collapse = " + "))),
  family = poisson, data = bat_data
)
summary(m_rich_pois)
check_overdispersion(m_rich_pois)


# ------------------------------------------------------------
# 3.2) Negative Binomial GLM
# ------------------------------------------------------------
m_rich_nb <- glm.nb(
  as.formula(paste(response_richness, "~", paste(predictors, collapse = " + "))),
  data = bat_data
)
summary(m_rich_nb)
check_overdispersion(m_rich_nb)


# ------------------------------------------------------------
# 3.3) GLMM (if random effects are available)
# ------------------------------------------------------------
if (random_site %in% colnames(bat_data)) {
  m_rich_glmm <- glmmTMB(
    as.formula(paste(response_richness, "~",
                     paste(predictors, collapse = " + "),
                     "+ (1 |", random_site, ")")),
    family = nbinom2, data = bat_data
  )
  summary(m_rich_glmm)
  AIC(m_rich_pois, m_rich_nb, m_rich_glmm)
} else {
  AIC(m_rich_pois, m_rich_nb)
}


# ------------------------------------------------------------
# 3.4) Diagnostics
# ------------------------------------------------------------
par(mfrow = c(2,2))
plot(m_rich_nb)
which(cooks.distance(m_rich_nb) > 4 / nrow(bat_data))

res.model <- simulateResiduals(m_rich_glmm)
plot(res.model)
testDispersion(res.model)
testZeroInflation(res.model)
testResiduals(res.model)


# ------------------------------------------------------------
# 3.5) Visualization
# ------------------------------------------------------------
plot_model(m_rich_nb, type = "est", transform = NULL, show.values = TRUE, show.p = TRUE) +
  theme_bw() + labs(title = "Species richness model", y = "Predictors", x = "Estimates")



# ==========================================================
# 04) EVENNESS / DIVERSITY MODELING
# ==========================================================


# ------------------------------------------------------------
# 4.1) Gamma GLM
# ------------------------------------------------------------
m_even <- glm(
  as.formula(paste(response_evenness, "~", paste(predictors, collapse = " + "))),
  family = Gamma(link = "log"), data = bat_data)
summary(m_even)


# ------------------------------------------------------------
# 4.2) Visualization
# ------------------------------------------------------------
plot_model(m_even, type = "est", transform = NULL, show.values = TRUE, show.p = TRUE) +
  theme_bw() + labs(title = "Evenness / diversity model", y = "Predictors", x = "Estimates")


# ==========================================================
# 05) NON-LINEAR EFFECTS (GAM)
# ==========================================================


# ------------------------------------------------------------
# 5.1) GAM for activity
# ------------------------------------------------------------
m_act_gam <- gam(
  as.formula(paste(response_activity, "~",
                   paste(paste0("s(", predictors, ")"), collapse = " + "))),
  family = nb(), data = bat_data
)
summary(m_act_gam)


# ------------------------------------------------------------
# 5.2) Model comparison
# ------------------------------------------------------------
AIC(m_act_nb, m_act_gam)


# ------------------------------------------------------------
# 5.3) Smooth term visualization
# ------------------------------------------------------------
par(mfrow = c(1, length(predictors)))
plot(m_act_gam, shade = TRUE)



# ==========================================================
# 06) SPECIES-LEVEL ANALYSIS (LONG FORMAT + GLMM)
# ==========================================================


# ------------------------------------------------------------
# 6.1) Reshape to long format
# ------------------------------------------------------------
data_long <- bat_data %>%
  rowwise() %>%
  mutate(across(all_of(species_columns), as.numeric)) %>%
  pivot_longer(cols = all_of(species_columns),
               names_to = "species",
               values_to = "species_activity") %>%
  mutate(non_activity = !!sym(response_activity) - species_activity) %>%
  ungroup()


# ------------------------------------------------------------
# 6.2) Multi-species GLMM (binomial)
# ------------------------------------------------------------
m_multi <- glmmTMB(as.formula(paste("cbind(species_activity, non_activity) ~", 
                   paste(predictors, collapse = " + "), "+ (1 | species)")),
  family = binomial, data = data_long)
summary(m_multi)


# ------------------------------------------------------------
# 6.3) Diagnostics
# ------------------------------------------------------------
res.model <- simulateResiduals(m_multi)
plot(res.model)
testDispersion(res.model)
testZeroInflation(res.model)
testResiduals(res.model)


# ------------------------------------------------------------
# 6.4) Visualization
# ------------------------------------------------------------
plot_model(m_multi, type = "est", transform = NULL, show.values = TRUE, show.p = TRUE) +
  theme_bw() + labs(title = "Multi-species response model", y = "Predictors", x = "Estimates")


# ==========================================================
# 07) FUNCTIONAL GUILD ANALYSIS (OPTIONAL)
# ==========================================================

# Example structure (adapt to your dataset):
# data$guild_1_activity <- data$SPECIES_A + data$SPECIES_B
# data$guild_2_activity <- data$SPECIES_C + data$SPECIES_D

# m_guild1 <- glm.nb(guild_1_activity ~ predictor1 + predictor2, data = data)
# summary(m_guild1)


# ==========================================================
# 08) SAVE MODEL OUTPUTS
# ==========================================================

saveRDS(list(
  activity_nb  = m_act_nb,
  richness_nb  = m_rich_nb,
  evenness_glm = m_even,
  activity_gam = m_act_gam,
  multi_species = m_multi),
  file.path(output_dir, "models_03_univariate.rds")
)

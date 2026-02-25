# ============================================================
# Author: Alberto Israel
# 00_setup.R
# Purpose: Install and load all required packages
# Last update: 2026-02-24
# ============================================================

# List of required packages
packages <- c(
  "tidyverse", "vegan", "ade4", "factoextra",
  "MASS", "lme4", "glmmTMB", "DHARMa", "performance",
  "visreg", "sjPlot", "MuMIn", "mgcv", "corrplot",
  "indicspecies", "car", "cooccur", "betapart", "igraph",
  "here", "RColorBrewer")

# Install missing packages
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))

# Session info for reproducibility
cat("R version:", as.character(getRversion()), "\n")
cat("Packages loaded successfully.\n")
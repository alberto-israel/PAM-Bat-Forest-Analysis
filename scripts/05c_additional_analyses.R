# ============================================================
# 05_additional_analyses.R
# Author: Alberto Israel
# Purpose: Additional fine-scale community analyses
# Section c: Network analysis of bat community based on site x species matrix
# Notes:
#   - Designed to integrate with 01_data_import.R
#   - User must define site column, species columns
# Last update: 2026-02-24
# ============================================================

# ------------------------------------------------------------
# 0) Setup project environment
# ------------------------------------------------------------
source(here::here("scripts", "00_setup.R"))

data_dir   <- file.path(here::here(), "data")
output_dir <- file.path(here::here(), "outputs")

library(igraph)  # building and analyzing networks (nodes, edges, metrics)
library(dplyr)   # data manipulation
library(tidyr)   # data reshaping
library(ggplot2)  # plotting
library(RColorBrewer) # color palettes for plots

# ------------------------------------------------------------
# Load prepared dataset
bat_data <- read_csv(file.path(output_dir, "bat_data_cleaned.csv"))
# ------------------------------------------------------------



# ================================
# 1) User inputs
# ================================
site_col <- "site"                  # Name of the column with site IDs
species_cols <- 4:19                # Columns containing species presence/absence or abundance
threshold_abundance <- 0            # Minimum value to consider presence (0 = presence/absence)
weighted <- FALSE                   # TRUE if you want edge weights based on abundance
# Optional: environmental grouping variable
env_group <- "habitat_type"         # e.g., Forest / Open, can be NA

#Pre-analysis column check
required_cols <- c(site_col, env_group)
missing_cols <- setdiff(required_cols, colnames(bat_data))
if(length(missing_cols) > 0) warning(paste("Missing columns:", paste(missing_cols, collapse=", ")))

if(any(species_cols > ncol(bat_data))) stop("species_cols exceed dataset columns")

##Pre-analysis variable check
if(is.na(env_group)) {
  env_group <- NULL
}

# ================================
# 2) Prepare site x species matrix
# ================================
mat <- bat_data[, species_cols]
rownames(mat) <- bat_data[[site_col]]

# Convert to presence/absence
mat_pa <- ifelse(mat > threshold_abundance, 1, 0)

# Check if matrix is empty
if(nrow(mat_pa) == 0 | ncol(mat_pa) == 0) stop("Species matrix is empty after filtering")


# ================================
# 3) Build network
# ================================
# Compute species co-occurrence matrix (simple correlation)
co_matrix <- if(weighted){
  cor(mat, use="pairwise.complete.obs", method="pearson")
} else {
  cor(mat_pa, use="pairwise.complete.obs", method="pearson")
}

# Remove self-links
diag(co_matrix) <- 0

# Check if all zeros
if(all(co_matrix == 0)) warning("All correlations are zero: network will be empty")

# Convert to igraph object
g <- graph_from_adjacency_matrix(co_matrix, mode = "undirected", weighted = weighted, diag = FALSE)


# ================================
# 4) Network metrics
# ================================
network_metrics <- data.frame(
  Degree = degree(g),
  Strength = if(weighted) strength(g, weights = E(g)$weight) else degree(g),
  Betweenness = betweenness(g, weights = if(weighted) 1/E(g)$weight else NULL),
  Closeness = closeness(g)
)

print(network_metrics)

# ================================
# 5) Plot network
# ================================
if(vcount(g) > 0){
  deg <- degree(g)
  cols <- colorRampPalette(brewer.pal(9,"YlGnBu"))(max(deg)+1)
  V(g)$color <- cols[deg+1]
  V(g)$size <- 5 + deg*2
  V(g)$label <- V(g)$name
  
  plot(g, layout = layout_with_fr, edge.width = if(weighted) E(g)$weight*2 else 1,
       vertex.label.cex = 0.8, vertex.label.color = "black",
       main = "Bat species network")
} else {
  warning("Network is empty: no plot generated")
}

# ================================
# 6) Optional: split by environmental group
# ================================
if(!is.null(env_group)){
  groups <- na.omit(unique(bat_data[[env_group]]))
  for(gp in groups){
    mat_gp <- mat_pa[bat_data[[env_group]] == gp, ]
    if(nrow(mat_gp) == 0 | ncol(mat_gp) == 0) next
    co_matrix_gp <- cor(mat_gp, use="pairwise.complete.obs")
    diag(co_matrix_gp) <- 0
    g_gp <- graph_from_adjacency_matrix(co_matrix_gp, mode="undirected", diag=FALSE)
    
    png(file.path(output_dir, paste0("network_", gp, ".png")), width=800, height=600)
    if(vcount(g_gp) > 0){
      plot(g_gp, layout = layout_with_fr, main = paste("Network:", gp))
    }
    dev.off()
  }
}

# ================================
# 7) Save results
# ================================
# 7.1) Save network metrics
write_csv(network_metrics, file.path(output_dir, "network_metrics.csv"))
message("Network metrics saved to outputs")

# 7.2) Save network plot
png(file.path(output_dir, "network_overall.png"), width=800, height=600)
if(vcount(g) > 0){
  plot(g, layout = layout_with_fr, edge.width = if(weighted) E(g)$weight*2 else 1,
       vertex.label.cex = 0.8, vertex.label.color = "black",
       main = "Bat species network")
}
dev.off()


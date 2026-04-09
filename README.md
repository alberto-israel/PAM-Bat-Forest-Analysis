#PAM Bat Forest Analysis

## Description
This repository contains a set of R scripts for analyzing bat ecology using passive acoustic monitoring data at the site × species level. The workflow focuses on acoustic activity (total and species-specific), species richness, diversity (Simpson/Shannon), and community composition. The scripts are designed for flexible integration with user-provided datasets and include both univariate and multivariate approaches, as well as additional analyses to explore finer community patterns.

### Scripts overview
1. **01_data_import.R**  
   Handles dataset import, cleaning, and initial formatting. Prepares a standardized `bat_data` object for downstream analyses.

2. **02_multivariate_exploration.R**  
   Performs exploratory multivariate analyses of community composition, including PCA and NMDS with ordination plots. Allows visual inspection of patterns before formal modeling. Allows visual inspection of patterns before formal modeling.

3. **03_modeling_univariate.R**  
   Implements GLM, GLMM, and GAM approaches to test ecological hypotheses on bat activity, species richness, and evenness as functions of environmental predictors. Generates coefficient plots and response visualizations.

4. **04_multivariate_inference.R**  
   Conducts inferential analyses at the community level, including PERMANOVA, RDA/CCA, variance partitioning, and indicator species analysis, to detect which environmental variables drive community composition and which species are characteristic of particular habitats

5. **05_additional_analyses.R**  
   Contains complementary analyses that reveal finer community responses, including:  
   - **05a_cooccurrence**: evaluates pairwise species co-occurrence patterns in different habitat types.  
   - **05b_beta_diversity**: computes site-to-site dissimilarities and partitions beta diversity into turnover and nestedness components.  
   - **05c_network_analysis**: builds species interaction networks to visualize potential associations and modularity within the community.

6. **06_sampling_completeness.R**
     Evaluates sampling completeness and adequacy using: species accumulation curves (spatial and temporal), non-parametric richness estimators (Chao1, Jackknife), coverage-based rarefaction (iNEXT), stratified completeness analysis by habitat type.

## Folder structure
- `data/`: raw and example datasets (e.g., `bat_data.csv`)  
- `scripts/`: R scripts (01-05)  
- `outputs/`: graphs, tables, and result files  
- `README.md`: this file

## Requirements
- R ≥ 4.0.0
- RStudio (recommended)
- Required packages are automatically installed via `00_setup.R`
See `scripts/00_setup.R` for the complete list of dependencies.

## How to run
This project uses an `.Rproj` file together with the `here` package

To run the analysis correctly:
1. Download or clone the repository.
2. Open the project by double-clicking `PAM_Bat_Forest_Analysis.Rproj`.
3. Open RStudio (the working directory will be set automatically).
4. Place your dataset in `data/bat_data.csv`.
5. Run the scripts in sequential order: 01 → 05.
6. Results will be saved in the `outputs/` folder.

## Citation & Context
This workflow was developed as part of my MSc thesis in Nature Sciences 
at Sapienza University of Rome (2025). 
The pipeline was designed to be generalizable and reproducible, allowing 
users to adapt it to their own passive acoustic monitoring datasets.

If you use or adapt this code, please cite this repository.


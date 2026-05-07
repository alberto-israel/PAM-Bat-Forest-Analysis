# PAM Bat Forest Analysis

## Description

This repository implements a reproducible R workflow to analyse bat communities from passive acoustic monitoring (PAM) data. The framework translates acoustic detections into ecological inference on how spatial covariates shape bat activity, species diversity, and community composition across scales.
The pipeline supports both raw Kaleidoscope Pro outputs and pre-aggregated site × species matrices, ensuring flexibility across datasets while maintaining a consistent analytical structure. It integrates data preparation, statistical modelling, multivariate community analysis, sampling completeness assessment, and environmental covariate extraction within a unified workflow.
The workflow is designed as a generalizable template for structuring ecological analyses of acoustic data, from raw detections to community-level inference.

---

## Overview

The pipeline is organized in three sections:

**Section 0 — Setup & Data Preparation** (run once before analysis)

- **`00a_setup.R`** — Installs and loads all required R packages.

- **`00b_data_preparation.R`** *(optional — only if starting from raw Kaleidoscope output)* — Imports raw Kaleidoscope `id.csv` files from multiple sites, performs quality control, assigns nocturnal sessions, and produces three analysis-ready matrices: site × species (aggregated), site × day × species, and site × hour × day × species. Also extracts site coordinates from Kaleidoscope summary files. All outputs are saved to `outputs/`.

**Section 1 — Site × Species Analysis** (scripts 01–06)

- **`01_data_import.R`** — Imports the site × species matrix, cleans data, and computes community-level response variables (total activity, guild-specific activity, species richness, Simpson and Shannon diversity indices).

- **`02_multivariate_exploration.R`** — Exploratory multivariate analyses of community composition, including PCA and NMDS with ordination plots. Allows visual inspection of patterns before formal modeling.

- **`03_modeling_univariate.R`** — GLM, GLMM, and GAM approaches to test ecological hypotheses on bat activity, species richness, and evenness as functions of environmental predictors. Generates coefficient plots and response visualizations.

- **`04_multivariate_inference.R`** — Inferential analyses at the community level: PERMANOVA, RDA/CCA, variance partitioning, and indicator species analysis. Identifies which environmental variables drive community composition and which species are characteristic of particular habitats.

- **`05_additional_analyses.R`** — Complementary analyses revealing finer community responses:
  - **05a** — Pairwise species co-occurrence patterns across habitat types.
  - **05b** — Site-to-site dissimilarity with beta diversity partitioned into turnover and nestedness components.
  - **05c** — Species interaction networks to visualize potential associations and community modularity.

- **`06_sampling_completeness.R`** — Sampling completeness and adequacy assessment: species accumulation curves (spatial and temporal), non-parametric richness estimators (Chao1, Jackknife), coverage-based rarefaction (iNEXT), stratified completeness analysis by habitat type.

**Section 2 — Environmental Covariates** (scripts 07–08)

- **`07_temperature.R`** — Extracts nocturnal temperature data from Kaleidoscope summary files (`S*.txt`). Produces hourly, daily, and site-level temperature tables and joins them with the detection matrices.

- **`08_moon.R`** — Computes lunar covariates (moon altitude, azimuth, illuminated fraction, phase, and brightness) for each sampled site × night × hour combination using `suncalc`. Joins results with the hourly detection matrix. Moon brightness (illuminated fraction × visibility above horizon) is the primary ecologically relevant output.

---

## Folder structure

```
PAM-Bat-Forest-Analysis/
├── data/
│   ├── example/
│   │   └── bat_data_site_species.csv   ← example site × species matrix (use for scripts 01–06)
│   └── raw_pam_outputs/                ← place your Kaleidoscope output here (see below)
│       ├── Output_Punto_1/
│       │   └── id.csv
│       ├── Output_Punto_2/
│       │   └── id.csv
│       ├── ...
│       ├── S1.txt
│       ├── S2.txt
│       └── ...
├── scripts/
│   ├── 00a_setup.R
│   ├── 00b_data_preparation.R
│   ├── 01_data_import.R
│   ├── 02_multivariate_exploration.R
│   ├── 03_modeling_univariate.R
│   ├── 04_multivariate_inference.R
│   ├── 05a_additional_analyses.R
│   ├── 05b_additional_analyses.R
│   ├── 05c_additional_analyses.R
│   ├── 06_sampling_completeness.R
│   ├── 07_temperature.R
│   └── 08_moon.R
├── outputs/                            ← generated automatically; all results saved here
├── README.md
├── LICENSE
└── PAM_Bat_Forest_Analysis.Rproj
```

### How to prepare `raw_pam_outputs/`

Kaleidoscope produces one output folder per recording site. Each folder must contain an `id.csv` file. The summary `.txt` files (one per site) must be placed directly inside `raw_pam_outputs/`. The expected naming convention is:

```
raw_pam_outputs/
├── Output_Punto_1/
│   └── id.csv
├── Output_Punto_2/
│   └── id.csv
├── S1.txt
└── S2.txt
```

Folder and file names can differ — `00b` extracts the site ID automatically from the first number in the folder name. If extraction fails, it falls back to the folder name itself.

---

## Requirements

- R ≥ 4.0.0
- RStudio (recommended)
- All required packages are installed automatically via `00a_setup.R`. See `scripts/00a_setup.R` for the full list of dependencies.

---

## How to run

This project uses an `.Rproj` file together with the `here` package for portable path management.

**Option A — Starting from raw Kaleidoscope output**

1. Clone or download the repository.
2. Open `PAM_Bat_Forest_Analysis.Rproj` in RStudio.
3. Place your Kaleidoscope output folders and summary `.txt` files in `data/raw_pam_outputs/` (see structure above).
4. Run `00a_setup.R` to install and load packages.
5. Run `00b_data_preparation.R` to generate analysis matrices and extract site coordinates (saved to `outputs/`).
6. Run `07_temperature.R` and `08_moon.R` to add environmental covariates.
7. Run scripts `01` → `06` in order. Results are saved to `outputs/`.

**Option B — Starting from a pre-aggregated site × species matrix**

1. Clone or download the repository.
2. Open `PAM_Bat_Forest_Analysis.Rproj` in RStudio.
3. Place your matrix in `data/example/bat_data_site_species.csv` (see the example file for the expected format).
4. Run `00a_setup.R` to install and load packages.
5. Run scripts `01` → `06` in order. Results are saved to `outputs/`.

---

## Citation & Context

This workflow was developed as part of my MSc thesis in Nature Sciences at Sapienza University of Rome (2025).
If you use or adapt this code, please cite this repository.

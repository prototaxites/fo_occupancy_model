# Future Oak - Occupancy Model and Amplicon Sequencing Analysis

This repository contains the raw data and analysis code for the paper "Environment and disease have tissue-specific effects on the tree microbiome" (Downie et al., *Environment and disease have tissue-specific effects on the tree microbiome*, bioRxiv, 2025.08.15.670310, https://doi.org/10.1101/2025.08.15.670310).

This project uses a `{targets}` pipeline to process raw amplicon sequencing data and builds a three-level occupancy model describing oak microbiota distribution across Britain. The occupancy model is implemented in Stan.

## System Requirements

### Operating Systems
The software has been tested on:
- **Linux** (Ubuntu 22.04)

### Software Dependencies

The analysis requires:

- **R 4.3.3** or later
- **Quarto 1.5.36** (for rendering reports)
- **cmdstan 2.34.1** (installed via R package `cmdstanr`)

#### R Packages (Core)
- `targets` (workflow management)
- `tarchetypes` (target recipes)
- `tidyverse` (data manipulation and visualization)
- `phyloseq` (microbiome data analysis)
- `brms` (Bayesian regression)
- `cmdstanr` (Stan interface, R ≥ 4.0)
- `posterior` (posterior distribution handling)
- `future` and `future.batchtools` (parallel computing)

#### R Packages (Analysis)
- `microViz`, `microbiome`, `speedyseq` (microbiome analysis)
- `vegan` (community ecology)
- `sf`, `sp` (spatial data)
- `ggdist`, `ggVennDiagram`, `patchwork`, `ComplexHeatmap` (visualization)
- `metacoder` (taxonomy visualization)
- `quarto` (rendering)
- `janitor`, `here` (utilities)

### Hardware Requirements
- **Minimum**: 8 GB RAM for basic analysis
- **Recommended**: 200+ GB RAM for full occupancy model fitting (3-level hierarchical Bayesian models across 6 tissue/marker combinations)
- **Optimal**: High-performance computing (HPC) cluster with SLURM scheduler for parallel model runs

## Installation Guide

### Option 1: Local Installation (Recommended for Development)

**Prerequisites**: Install [R 4.3.3+](https://cran.r-project.org/) and [Quarto](https://quarto.org/docs/get-started/)

**Steps**:

1. Clone the repository:
   ```bash
   git clone https://github.com/prototaxites/fo_occupancy_model.git
   cd fo_occupancy_model
   ```

2. Install R dependencies:
   ```r
   # Install pak for faster dependency resolution
   install.packages("pak", repos = sprintf(
     "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
     .Platform$pkgType, R.Version()$os, R.Version()$arch
   ))
   
   # Add custom repositories
   pak::repo_add(Stan = "https://mc-stan.org/r-packages/",
                 microViz = "https://david-barnett.r-universe.dev")
   
   # Install all dependencies
   pak::pak(c('targets', 'tarchetypes', 'quarto', 'qs',
              'tidyverse', 'here', 'janitor', 
              'microViz', 'metacoder', 'vegan',
              'posterior', 'brms', 
              'sf', 'future', 'future.batchtools', 
              'ggdist', 'ggVennDiagram', 'patchwork',
              'RcppEigen', 'bayesm', 'BH'))
   
   pak::pkg_install(c('phyloseq', 'microbiome', 'microViz',
                      'cmdstanr', 'ComplexHeatmap', 'mikemc/speedyseq'))
   ```

3. Install cmdstan:
   ```r
   cmdstanr::install_cmdstan(version = "2.34.1")
   ```

### Option 2: Docker/Singularity (Recommended for Reproducibility)

Pre-built containers include all dependencies and cmdstan:

```bash
# Build image from Dockerfile
docker build -t fo_occupancy_model:latest -f container/Dockerfile .

# Run analysis in container
docker run -v $(pwd):/work -w /work \
  fo_occupancy_model:latest \
  Rscript -e "targets::tar_make()"
```

### Full Pipeline Run

**Expected time**: 48+ hours on an HPC cluster with 64 CPU cores and 100+ GB RAM

```r
library(targets)
tar_make()  # Run entire pipeline
```

**Expected Output**:
- 6 compiled Stan models (one per tissue/marker combination)
- Posterior draws and predictions for each model
- Summary statistics and diagnostics
- Landscape-scale ASV occupancy predictions
- Model convergence diagnostics (Rhat, ESS)

## Instructions for Use

### Running the Full Analysis

1. **Prepare your environment**:
   ```r
   library(targets)
   tar_source()  # Load custom functions
   ```

2. **View the pipeline**:
   ```r
   tar_visnetwork()  # Interactive DAG visualization
   tar_manifest()    # Table of all targets
   ```

3. **Inspect results**:
   ```r
   # Load a completed target
   tar_load(model_ITS_Leaf)
   model_ITS_Leaf$summary()  # Model diagnostics
   
   # Load predictions
   tar_load(landscape_predictions_psi_ITS_Leaf)
   ```

### Running Downstream Analyses

All downstream analyses are in Quarto markdown files:

```bash
# Render all reports
quarto render reports/

# Render individual reports
quarto render reports/combined_analysis.qmd
quarto render reports/model_analysis.qmd
quarto render reports/alpha_diversity.qmd
```

**Key output files**:
- `reports/combined_analysis.html`: All main figures and supplementary figures
- Stan model files: `Stan/sgcp_occupancy_envsplit.stan` (main model)
- Model diagnostics: Accessible via `tar_load(model_diagnostics_[MARKER]_[TISSUE])`

**Code repositories**:
- Amplicon analysis (this repo): https://github.com/prototaxites/fo_occupancy_model
- Functional analysis: https://github.com/prototaxites/fo_functional_analysis

## References

For detailed methodology, see:
- Downie et al. (2025). "Environment and disease have tissue-specific effects on the tree microbiome." *bioRxiv*. https://doi.org/10.1101/2025.08.15.670310
- Methods section covers: site selection, sampling, DNA extraction, amplicon sequencing, metagenome processing, occupancy modeling, and statistical analysis.

---

**Funding**: This work was funded by UK Research and Innovation's Strategic Priorities Fund (BBSRC grant BB/T01069X/1) under the FUTURE OAK project.

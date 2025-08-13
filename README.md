# Future Oak - Occupancy Model and amplicon sequencing analysis

This repository contains the raw data and analysis code for the paper "Environment and disease have tissue-specific effects on the tree microbiome" (Downie et al., in prep).

The repository contains a `{targets}` pipeline, which processes the raw amplicon sequencing data in `data/`, 
and builds the occupancy model as described in the paper. This pipeline can be run using the folling R code:

```
library(targets)
targets::tar_make()
```

The occupancy model itself is written in Stan, and is available in the file `Stan/sgcp_occupancy_envsplit.stan`.

All downstream analyses used in the paper are coded as Quarto markdown files, and are available in the `reports` directory. The primary analyses are all
available in the `reports/combined_analyses.qmd` file.

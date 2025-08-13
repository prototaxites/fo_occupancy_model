#!/usr/bin/env Rscript

install.packages("pak", repos = sprintf(
  "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
  .Platform$pkgType,
  R.Version()$os,
  R.Version()$arch
))

pak::repo_add(.list = c(Stan = "https://mc-stan.org/r-packages/", 
                microViz = "https://david-barnett.r-universe.dev")
             )
			 
if(R.Version()$arch == "aarch64") {
	pak::repo_add(.list = c(RHUB = "https://raw.githubusercontent.com/r-hub/repos/main/ubuntu-22.04-aarch64/4.3"))
}

install.packages(c("rncl", "GA"))

pak::pak(c('targets', 'tarchetypes',  'quarto', 'qs',
                   'tidyverse', 'here', 'janitor', 
                   'microViz', 'metacoder', 'vegan',
                   'posterior', 'brms', 
                   'sf', 'future', 'future.batchtools', 
                   'ggdist', 'ggVennDiagram','patchwork',
                   'RcppEigen', 'bayesm', 'BH'
                   ))

pak::pkg_install(c('phyloseq','microbiome', 'microViz', 
                   'cmdstanr', 'ComplexHeatmap', 'mikemc/speedyseq'
))
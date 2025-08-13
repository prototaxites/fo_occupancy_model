#!/usr/bin/env Rscript

## Builds the Quarto documents in the project
## I tried having the {targets} pipeline do this, but some combination of 
## tar_make_future + Singularity + Slurm makes the jobs hang!

library(tidyverse)
library(glue)
library(quarto)

wd <- here::here()

parameterised_build <- function(quarto_file, prefix, execute_params, output_dir) {
  marker <- execute_params$marker
  tissue <- execute_params$tissue
  output_file <- glue::glue("{prefix}_{marker}_{tissue}.html")

  quarto_render(
    input = quarto_file,
    output_format = "html",
    output_file = output_file,
    execute_params = execute_params,
    execute_dir = wd
  )

  file.rename(output_file, glue::glue("{output_file}"))
}

map_df <- expand_grid(marker = c("ITS", "16S"), 
                      tissue = c("Leaf", "Stem", "Soil")
)
map_list <- split(map_df, 
                  f = seq_len(nrow(map_df))
                  ) |> 
            map(\(x) list(marker = x$marker, tissue = x$tissue))

## Build diagnostic docs
walk(map_list, \(x) parameterised_build("model_diagnostics.qmd", "diagnostics", x, "reports"))

## Build analysis docs
walk(map_list, \(x) parameterised_build("model_analysis.qmd", "analysis", x, "reports"))

## Build alpha diversity
quarto_render(
  input = "alpha_diversity.qmd",
  output_format = "html",
  output_file = "alpha_diversity.html",
  execute_dir = wd
)
## Build beta diversity

quarto_render(
  input = "beta_diversity.qmd",
  output_format = "html",
  output_file = "beta_diversity.html",
  execute_dir = wd
)


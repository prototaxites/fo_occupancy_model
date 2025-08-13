library(targets)
library(tarchetypes)
library(future)
library(future.batchtools)
library(tidyverse)

## Load functions
tar_source()

# Set target-specific options such as packages.
tar_option_set(packages = c("tidyverse", "phyloseq",  "brms", "cmdstanr", "posterior", "microViz"),
               format = "qs",
               memory = "transient",
               garbage_collection = TRUE,
               storage = "worker",
               retrieval = "worker",
               error = "abridge"
)


plan(batchtools_slurm, 
     template = "tmpl/future.tmpl", 
     resources = list(ncpus = 1, walltime = "5-00:00:0", memory = "100G", partition = "large"),
     workers = 10
)

model_resources <- tar_resources(
  future = tar_resources_future(
    plan = tweak(
      batchtools_slurm,
      template = "tmpl/future_overlay.tmpl",
      resources = list(ncpus = 64, memory = "1028900", walltime = "5-00:00:0", partition = "large")
    )
  )
)

model_subsets <- tidyr::expand_grid(
  marker = c("ITS", "16S"),
  tissue = c("Leaf", "Stem", "Soil")) |> 
  mutate(phyloseq = ifelse(marker == "ITS", 
                           here::here("data/its_microbiome_phyloseq.rds"), 
                           here::here("data/16s_microbiome_phyloseq.rds")),
         filter_thresh = ifelse(marker == "ITS", 100, 50)
  )

list(
  ## Preprocesses tree and site metadata
  tar_file(name = stan_model, 
           command = "Stan/sgcp_occupancy_envsplit.stan"),
  tar_file(name = stan_postpred, 
           command = "Stan/sgcp_occupancy_envsplit_postpred.stan"),
  tar_file_read(name = site_df,
                command = here::here("data/site_metadata.csv"),
                read = read_csv(!!.x, col_types = cols()) |> arrange(Site)
  ),
  tar_file_read(name = tree_df,
                command = here::here("data/trees.csv"),
                read = read_csv(!!.x, col_types = cols())
  ),
  tar_file_read(name = decline_df,
                command = here::here("data/decline.csv"),
                read = read_csv(!!.x, col_types = cols()) |>
                  arrange(Site, Tree) |>
                  mutate(PDI_rescaled = scale(PDI_unscaled)[,1])
  ),
  tar_file_read(name = pheno_df,
                command = here::here("data/phenotypes.csv"),
                read = read_csv(!!.x, col_types = cols()) |>
                  janitor::clean_names(case = "upper_camel") |>
                  arrange(Site, Tree)
  ),
  tar_file_read(name = landscape_df, 
                command = here::here("data/uk_wide_env/uk_wide_env.shp"),
                read = read_landscape(!!.x,
                                      vars = c("groundfrost", "hurs", "pv", "rainfall", 
                                               "sfcWind", "sun", "tas", "Ntot", "SOx", "CaMg"
                                               )
                                      ) |> 
                  filter(Northing <= 710000) |> 
                  filter(!(Easting < 2e5 & between(Northing, 450000, 600000)))
  ),
  tar_target(name = landscape_sample_df,
             command = data.frame(Disease = rep(0, nrow(landscape_df)))
  ),
  tar_target(name = env_vars, 
             command = c("groundfrost", "hurs", "pv", "rainfall", 
                         "sfcWind", "sun", "tas", "Ntot", "SOx", "CaMg")
  ),
  tar_target(name = sample_vars, 
             command = "Disease"
  ),
  tar_map(
    ## build dataframe of tissues/kingdom to map over
    values = model_subsets,
    names = c("marker", "tissue"),
    tar_file_read(name = microbiome, 
                  command = phyloseq, 
                  read = readRDS(!!.x)
    ),
    tar_target(name = microbiome_filtered, 
               command = microbiome |>
                 subset_taxa(Kingdom %in% c("Fungi", "Bacteria")) |>
                 subset_taxa(Order != "Chloroplast") |>
                 subset_taxa(Family != "Mitochondria") |>
                 ps_filter(Tissue == tissue) |>
                 ps_filter(!grepl("E", Tree)) |>
                 filter_taxa(\(x) any(x > filter_thresh), prune = TRUE) |> 
                 {\(x) prune_samples(sample_sums(x) != 0, x)}()
    ),
    tar_target(name = asv_table,
               command = microbiome_filtered |> 
                 tax_table() |>
                 as.data.frame() |> 
                 rownames_to_column("ASV_ID") |> 
                 as_tibble() |> 
                 mutate(ASV_numeric = as.numeric(as.factor(ASV_ID))) |> 
                 arrange(ASV_numeric)
    ),
    tar_target(name = taxa,
               command = 1:ntaxa(microbiome_filtered)
    ),
    tar_target(name = microbiome_tj_corrected,
               command = ps_tagjump_correct(microbiome_filtered),
               resources = model_resources
    ),
    tar_target(name = stan_data,
               command = make_standata_envsplit(microbiome_tj_corrected,
                                                site_df,
                                                env_vars = env_vars,
                                                sample_vars = sample_vars
               )
    ),
    tar_target(name = model,
               command = cmdstan_fit_model(stan_file = stan_model,
                                           data = stan_data,
                                           chains = 4,
                                           iter_warmup = 1000,
                                           iter_sampling = 1000,
                                           parallel_chains = 4,
                                           threads_per_chain = 16,
                                           sig_figs = 9,
                                           cpp_options = list(stan_threads = TRUE),
                                           adapt_delta = 0.9,
                                           step_size = 0.05,
                                           init = 0.1,
                                           seed = 12345
               ), 
               resources = model_resources),
    tar_target(name = model_draws, 
               command = model$draws()
    ),
    tar_target(name = model_summary, 
               command = model$summary()
    ),
    tar_target(name = model_diagnostics, 
               command = model$diagnostic_summary()
    ),
    tar_target(name = model_postpreds,
               command = cmdstan_posterior_predict(
                 stan_file = stan_postpred,
                 data = stan_data,
                 draws = model_draws,
                 output_dir = NULL,
                 output_basename = NULL,
                 seed = 12345
                 ),
               resources = model_resources
    ),
    tar_target(name = model_postpred_richness,
               command = subset_draws(model_postpreds, 
                 c("Z_cond_richness", 
                   "U_cond_richness")
                )
    ),
    tar_target(name = loo_diag,
               command = model$loo()
    ),
    tar_target(name = landscape_predictions_psi,
               command = landscape_predictions(
                 draws = model_draws,
                 ref_df = site_df,
                 landscape_df = landscape_df,
                 sample_df = landscape_sample_df,
                 .fun = posterior_predict_phi,
                 chunk_size = 20
               ) |> reduce(cbind)
    )
  ),
  tar_target(
    name = phyloseq_its,
    command = phyloseq::merge_phyloseq(
      microbiome_tj_corrected_ITS_Leaf,
      microbiome_tj_corrected_ITS_Stem,
      microbiome_tj_corrected_ITS_Soil)
  ),
  tar_target(
    name = phyloseq_16s,
    command = phyloseq::merge_phyloseq(
      microbiome_tj_corrected_16S_Leaf,
      microbiome_tj_corrected_16S_Stem,
      microbiome_tj_corrected_16S_Soil)
  )
)

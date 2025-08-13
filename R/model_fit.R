cmdstan_fit_model <- function(stan_file,
                              data,
                              chains,
                              iter_warmup = 1000,
                              iter_sampling = 1000,
                              parallel_chains,
                              threads_per_chain = 1,
                              cpp_options = NULL,
                              sig_figs = 6,
                              inc_warmup = FALSE,
                              step_size = 1,
                              adapt_delta = 0.8,
                              init = 2,
                              seed = NULL) {
  
  if(!is.null(seed)) set.seed(seed)

  model <- cmdstanr::cmdstan_model(stan_file,
                                   cpp_options = cpp_options)
  fit <- model$sample(
    data = data,
    chains = chains,
    parallel_chains = parallel_chains,
    threads_per_chain = threads_per_chain,
    step_size = step_size,
    adapt_delta = adapt_delta,
    init = init,
    sig_figs = sig_figs
  )
  
  fit$draws(inc_warmup = inc_warmup)
  try(fit$sampler_diagnostics(inc_warmup = inc_warmup), silent = TRUE)
  try(fit$init(), silent = TRUE)
  try(fit$profiles(), silent = TRUE)
  fit
}

make_standata_multispp <- function(ps,
                                   site_data,
                                   env_vars,
                                   sample_vars,
                                   prior_sigma,
                                   prior_b,
                                   prior_corr) {

  N_species <- ntaxa(ps)
  N_replicates <- phyloseq::nsamples(ps)

  Y_df <- ps |> 
    psmelt() |> 
    ungroup() |> 
    select(OTU, Site, Tree, Abundance) |>
    arrange(OTU, Site, Tree) |>
    mutate(Site = as.numeric(as.factor(Site)),
           OTU = as.numeric(as.factor(OTU)),
           Rep = rep(1:N_replicates, times = N_species))
  
  Q <- Y_df |>
    group_by(Site, OTU) |>
    summarise(Q = ifelse(any(Abundance > 0), 1, 0)) |> 
    pivot_wider(values_from = Q, names_from = Site) |> 
    select(-OTU) |> as.matrix()
  
  Y <- select(Y_df, OTU, Rep, Abundance) |>
    pivot_wider(names_from = "Rep", values_from = "Abundance") |> 
    select(-OTU) |> as.matrix()
  
  site_df_scaled <- site_data |> 
    mutate(across(
      where(\(x) is.numeric(x)), 
      \(x) scale(x)[,1])
      )
  
  X_occ <- site_df_scaled |>
    mutate(Site = as.numeric(as.factor(Site))) |>
    arrange(Site) |>
    model.matrix(reformulate(env_vars), data = _)
  
  sample_data_processed <- sample_data(ps) |>
    as("data.frame") |> 
    rownames_to_column("SampleID") |>
    arrange(Site, Tree) |> 
    left_join(site_df_scaled)
  
  X_rep <- sample_data_processed |>
    model.matrix(c(env_vars, sample_vars), data = _)
  
  rep_group_sizes <- sample_data_processed |>
    arrange(Site) |>
    pull(Site) |>
    as.factor() |>
    as.numeric() |>
    table() |>
    as.vector()
  
  X_reads <- sample_data_processed |>
    model.matrix(c(env_vars, sample_vars), data = _)
  
  list(
    N_species = N_species,
    N_sites = nrow(X_occ),
    N_replicates = N_replicates,
    Y = Y,
    Q = Q,
    rep_group_sizes = rep_group_sizes,
    N_pred_occ = ncol(X_occ),
    X_occ = X_occ,
    N_pred_rep = ncol(X_rep),
    X_rep = X_rep,
    N_pred_reads = ncol(X_reads),
    X_reads = X_reads,
    grainsize = 1,
    prior_sigma = prior_sigma,
    prior_b = prior_b,
    prior_corr = prior_corr
  )
}

cmdstan_posterior_predict <- function(stan_file,
                              data,
                              draws,
                              output_dir = NULL,
                              output_basename = NULL,
                              seed = NULL) {
  
  if(!is.null(seed)) set.seed(seed)

  dummy_model <- cmdstan_model(stan_file, force_recompile = TRUE)
  gq <- dummy_model$generate_quantities(
    fitted_params = draws,
    data = data,
    output_dir = output_dir,
    output_basename = output_basename
  )
  gq$draws()
}

cmdstan_loo_moment_match <- function(stan_file,
                                     data,
                                     model_draws,
                                     cpus = 1) {
  
  dummy_model <- cmdstan_model(stan_file, force_recompile = TRUE)
  mod <- dummy_model$sample(data = data,
                            iter_warmup = 1,
                            iter_sampling = 1,
                            chains = chains, parallel_chains = parallel_chains)
  mod$.__enclos_env__$private$draws_ <- model_draws
  mod$loo(moment_match = TRUE, cores = cpus)
}

## ENVSPLIT
make_standata_envsplit <- function(ps,
                                   site_data,
                                   env_vars,
                                   sample_vars) {
  
  N_species <- ntaxa(ps)
  N_replicates <- phyloseq::nsamples(ps)
  
  Y_df <- ps |> 
    psmelt() |> 
    ungroup() |> 
    select(OTU, Site, Tree, Abundance) |>
    arrange(OTU, Site, Tree) |>
    mutate(Site = as.numeric(as.factor(Site)),
           OTU = as.numeric(as.factor(OTU)),
           Rep = rep(1:N_replicates, times = N_species))
  
  Q <- Y_df |>
    group_by(Site, OTU) |>
    summarise(Q = ifelse(any(Abundance > 0), 1, 0)) |> 
    pivot_wider(values_from = Q, names_from = Site) |> 
    select(-OTU) |> as.matrix()
  
  Y <- select(Y_df, OTU, Rep, Abundance) |>
    pivot_wider(names_from = "Rep", values_from = "Abundance") |> 
    select(-OTU) |> as.matrix()
  
  site_df_scaled <- site_data |> 
    mutate(across(
      where(\(x) is.numeric(x)), 
      \(x) scale(x)[,1])
    )
  
  X_env <- site_df_scaled |>
    mutate(Site = as.numeric(as.factor(Site))) |>
    arrange(Site) |>
    model.matrix(reformulate(env_vars), data = _)
  
  sample_data_processed <- sample_data(ps) |>
    as("data.frame") |> 
    rownames_to_column("SampleID") |>
    arrange(Site, Tree) |> 
    left_join(site_df_scaled)
  
  X_sample <- sample_data_processed |>
    model.matrix(reformulate(sample_vars, intercept = TRUE), data = _) |> 
    {\(x) x[,-1]}() |> 
    as.matrix()

  rep_group_sizes <- sample_data_processed |>
    arrange(Site) |>
    pull(Site) |>
    as.factor() |>
    as.numeric() |>
    table() |>
    as.vector()
  
  list(
    N_species = N_species,
    N_sites = nrow(X_env),
    N_replicates = N_replicates,
    Y = Y,
    Q = Q,
    rep_group_sizes = rep_group_sizes,
    N_pred_env = ncol(X_env),
    X_env = X_env,
    N_pred_sample = ncol(X_sample),
    X_sample = X_sample,
    grainsize = 1
  )
}

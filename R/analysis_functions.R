make_sigmas_df <- function(marker, tissue, level, env_vars, sample_vars) {
  if(level > 3) stop("Error: level must be between 0 and 3!")
  
  draws <- tar_read_raw(glue::glue("model_draws_{marker}_{tissue}"),
                        store = here::here("_targets")
  ) |>
    as_draws_rvars()
  
  if(level == 0) {
    sigmas_env <- subset_draws(draws, "sigmas_env")[[1]]
    sigmas_sample <- subset_draws(draws, "sigmas_sample")[[1]]
    mean <- NULL
    sigma <- c(sigmas_env, sigmas_sample)
    
    var_names <- c(env_vars, sample_vars)
  } else if(level > 0) {
    mean <- as.vector(subset_draws(draws, "b_env")[[1]][level,])
    sigma <- as.vector(subset_draws(draws, "sigmas_env_split")[[1]][,level])
    
    var_names <- c(env_vars)
    if(level > 1) {
      mean_sample <- as.vector(subset_draws(draws, "b_sample")[[1]][level - 1,])
      sigma_sample <- as.vector(subset_draws(draws, "sigmas_sample_split")[[1]][,level - 1])
      
      mean <- c(mean, mean_sample)
      sigma <- c(sigma, sigma_sample)
      
      var_names <- c(env_vars, sample_vars)
    }
  }
  
  tibble(marker = marker,
         tissue = tissue,
         var    = var_names,
         level  = level,
         mean   = mean,
         sigma  = sigma)
}

get_per_sample_raw_richness <- function(marker, tissue, split = TRUE, measures = NULL) {
  ps <- tar_read_raw(glue::glue("microbiome_tj_corrected_{marker}_{tissue}"))

  samplecounts <- phyloseq::sample_sums(ps)
  keep <- samplecounts > 0
  ps <- prune_samples(samples = keep, x = ps)
  sd <- sample_data(ps)
  
  estimate_richness(ps, split = split, measures = measures) |>
    mutate(Marker = marker, TotalReads = samplecounts[keep]) |> 
    cbind(sd)
}

model_richness_disease_differences <- function(marker, tissue, n_env_vars, chunk_size = NULL) {
  draws <- tar_read_raw(glue::glue("model_draws_{marker}_{tissue}"),
                        store = here::here("_targets")
  ) |> as_draws_rvars()
  
  env_data <- tibble(vars = 1:n_env_vars, value = 0) |> 
    pivot_wider(names_from = vars, values_from = value) |> 
    uncount(2) |> 
    model.matrix(~., data = _)
  
  sample_data <- matrix(c(0, 1), ncol = 1)
  
  ntaxa <- dim(draws$b_env_species)[3]
  taxa <- 1:ntaxa
  phi <- c(1, 1)

  if(is.null(chunk_size)) {
    taxa_chunks <- taxa
  } else {
    z <- seq_along(taxa)
    taxa_chunks <- split(taxa, ceiling(z/chunk_size))
    names(taxa_chunks) <- NULL
  }

  pred_richness <- map(taxa_chunks, \(x) posterior_predict_theta(x, draws, env_data, sample_data, phi)) |> 
    map(\(x) draws_of(x, with_chains = TRUE) |> 
          apply(c(1,2,3,4), \(y) rbinom(1, 1, y)) |> 
          rvar(with_chains = TRUE)
        ) |> 
    reduce(cbind) |>
    apply(1, rvar_sum)
  
  tibble(Marker = marker,
         Tissue = tissue,
         Disease = c(0, 1),
         pred = list_c(pred_richness))
}

make_richness_df <- function(marker, tissue, level) {
  if(level < 1) stop("Error: level must be between 1 and 3!")
  if(level > 2) stop("Error: level must be between 1 and 3!")
  draws <- tar_read_raw(glue::glue("model_postpred_richness_{marker}_{tissue}"),
                        store = here::here("_targets")
  ) |> 
    as_draws_rvars() 
  site_df <- tar_read(site_df)
  ps <- tar_read_raw(glue::glue("microbiome_tj_corrected_{marker}_{tissue}")) |> 
    ps_select(Site, Tree, Disease) |> 
    psmelt() |> 
    arrange(OTU, Site, Tree, Disease)
  
  if(level == 1) {
    ps <- summarise(ps, Abundance = sum(Abundance), .by = c(OTU, Site)) |> 
      summarise(Richness_obs = sum(Abundance > 0), .by = Site)
    richness_mod <- pluck(draws, "Z_cond_richness")
  } else {
    ps <- summarise(ps, Richness_obs = sum(Abundance > 0), .by = c(Site, Tree, Disease))
    richness_mod <- pluck(draws, "U_cond_richness")
  }
  
  ps <- mutate(ps, Join = as.character(1:n()))
  
  tibble(Richness_mod = richness_mod) |> 
    rownames_to_column("Join") |> 
    left_join(ps, by = "Join") |> 
    mutate(Richness_obs_diff = Richness_mod - Richness_obs,
           Marker = marker,
           Tissue = tissue) |> 
    left_join(site_df)
}

calculate_average_richness_differences <- function(marker, tissue) {
  richness_df <- make_richness_df(marker, tissue, level = 2)
  
  disease <- filter(richness_df, Disease == "DiseaseAOD") |> 
    pull(Richness_mod) |> as_draws_array() |> cbind() |> rvar()
  control <- filter(richness_df, Disease == "Control") |> 
    pull(Richness_mod) |> as_draws_array() |> cbind() |> rvar()
  
  tibble(
    Marker = marker,
    Tissue = tissue,
    Control = control,
    Disease = disease,
    Diff = Disease - Control,
    Model = "Occ_Avg"
  )
}

disease_intervals <- function(marker, tissue, chunk_size = NULL, n_env_vars = 10, tax_file = "data/asvs.taxonomy.sativa.tsv") {
  ps <- tar_read_raw(glue::glue("microbiome_tj_corrected_{marker}_{tissue}"),
                     store = here::here("_targets"))
  if(marker == "16S") ps <- fix_ps_tax(ps, tax_file)
  
  asv_table <- ps |> 
    tax_table() |> 
    as.data.frame() |> 
    rownames_to_column("ASV_ID") |> 
    as_tibble() |> 
    arrange(ASV_ID) |> 
    mutate(ASV_numeric = as.numeric(as.factor(ASV_ID)))
  
  model_draws <- tar_read_raw(glue::glue("model_draws_{marker}_{tissue}"),
                              store = here::here("_targets")) |> 
    as_draws_rvars()
  
  b_env_species <- pluck(model_draws, "b_env") + pluck(model_draws, "b_env_species")
  b_sample_species <- pluck(model_draws, "b_sample") + pluck(model_draws, "b_sample_species")
  
  n_pred_env <- dim(b_env_species)[2] - 1
  n_pred_sample <- dim(b_env_species)[2]
  n_taxa <- dim(b_env_species)[3]

  env_data <- tibble(vars = 1:n_env_vars, value = 0) |> 
    pivot_wider(names_from = vars, values_from = value) |> 
    uncount(2) |> 
    model.matrix(~., data = _)
  
  sample_data <- matrix(c(0, 1), ncol = 1)

  if(is.null(chunk_size)) {
    taxa_chunks <- 1:n_taxa
  } else {
    z <- seq_along(1:n_taxa)
    taxa_chunks <- split(z, ceiling(z/chunk_size))
    names(taxa_chunks) <- NULL
  }
  
  lambda <- map(taxa_chunks, \(x) posterior_predict_lambda(x, model_draws, env_data, sample_data, phi = c(1, 1))) |>
    reduce(cbind)
  
  t(lambda) |> 
    as.data.frame() |> 
    mutate(Marker = marker, Tissue = tissue) |> 
    rename(Control = `1`, DiseaseAOD = `2`) |> 
    rownames_to_column("ASV_numeric") |> 
    mutate(ASV_numeric = as.numeric(ASV_numeric)) |> 
    left_join(asv_table)
}

level_quantiles <- function(draws, level, env_vars) {
  n_env_vars <- length(env_vars)
  
  mean <- as.vector(subset_draws(draws, "b_env")[[1]][level, 2:(n_env_vars + 1)])
  sigma <- as.vector(subset_draws(draws, "sigmas_env_split")[[1]][2:(n_env_vars + 1),level])
  quantiles <- rvar_rng(rnorm, 10, 0, sigma) |> 
    quantile(probs = c(0.025, 0.975)) |> t()
  
  colnames(quantiles) <- c(
    glue::glue("lci"),
    glue::glue("uci")
  )
  
  as_tibble(quantiles) |> 
    mutate(coef = env_vars) |> 
    pivot_longer(cols = c("lci", "uci"), 
                 names_to = "ci",
                 values_to = "slope") |> 
    mutate(level = level)
}

expand_scaled_range <- function(site_data, var) {
  mean_var <- mean(pull(site_data, var))
  sd_var <- sd(pull(site_data, var))
  
  site_data <- mutate(site_data, across(where(is.numeric), \(x) scale(x)[,1]))
  
  expand_grid(coef = var, 
              min = min(pull(site_data, var)), 
              max = max(pull(site_data, var)),
              mean = mean_var,
              sd = sd_var,
              seq = seq(from = min, to = max, by = 0.1)) |> 
    select(-c(min, max)) |> 
    nest(.by = coef)
}

species_slopes_quantiles <- function(marker, tissue, level, env_vars, var_names, site_data){
  if(level < 1 | level > 3) stop("Error: level must be between 0 and 3!")
  
  draws <- tar_read_raw(glue::glue("model_draws_{marker}_{tissue}"),
                        store = here::here("_targets")
  ) |>
    as_draws_rvars()
  
  expansions <- map(env_vars, \(x) expand_scaled_range(site_data, x)) |> 
    list_rbind()
  
  res <- tibble(Tissue = tissue,
                Marker = marker,
                coef = env_vars,
                var_names = var_names) |> 
    left_join(expansions)
  
  quantiles <- map(1:level, \(x) level_quantiles(draws, x, env_vars)) |> 
    list_rbind() |> 
    left_join(res) |> 
    unnest(data)
  
  quantiles
}

representation_likelihood <- function(tax_level, sig_list, taxa_list, stan_file) {
  N_total <- nrow(taxa_list)
  N_pulled <- nrow(sig_list)
  
  counts_bag <- taxa_list |> 
    count(.data[[tax_level]], name = "n_bag")
  counts_sig <- sig_list |> 
    count(.data[[tax_level]], name = "n_pick")
  
  counts <- left_join(counts_bag, counts_sig) |> 
    mutate(n_pick = ifelse(is.na(n_pick), 0, n_pick))

  mod <- cmdstanr::cmdstan_model(stan_file)
  fit <- mod$sample(
    data = list(
      N = nrow(counts),
      y = counts$n_pick,
      mA = counts$n_bag,
      mB = N_total - counts$n_bag,
      n = N_pulled
    ),
    chains = 4,
    parallel_chains = 4
  )
  
  alpha <- fit$draws("alpha") |> 
    posterior::as_draws_rvars() |> 
    pluck("alpha")
  
  mutate(counts, alpha = alpha)
}

fix_ps_tax <- function(ps, tax_file) {
  old_order <- rownames(tax_table(ps))
  phyloplace_classifications <- read_tsv(tax_file, col_names = c("ASV_ID", "Lineage", "Support")) |> 
    mutate(Kingdom = str_extract(Lineage, "d__([[:alnum:]]+)\\;?", group = 1),
           Phylum = str_extract(Lineage, "p__([[:alnum:]_-]+)\\;?", group = 1),
           Class = str_extract(Lineage, "c__([[:alnum:]_-]+)\\;?", group = 1),
           Order = str_extract(Lineage, "o__([[:alnum:]_-]+)\\;?", group = 1),
           Family = str_extract(Lineage, "f__([[:alnum:]_-]+)\\;?", group = 1),
           Genus = str_extract(Lineage, "g__([[:alnum:]_-]+)\\;?", group = 1),
           Species = str_extract(Lineage, "s__([[:alnum:]_-]+)\\;?", group = 1)
    ) |> 
    select(ASV_ID, Kingdom:Species) |> 
    filter(ASV_ID %in% taxa_names(ps)) |> 
    column_to_rownames("ASV_ID") |> 
    as.matrix()
  
  phyloplace_classifications <- phyloplace_classifications[order(match(rownames(phyloplace_classifications), old_order)), , drop = FALSE]
  
  tax_table(ps) <- tax_table(phyloplace_classifications)
  
  return(ps)
}

estimate_colprods <- function(ps) {
  phylo_long_colprods <- ps |> 
    phyloseq::psmelt() |> 
    as_tibble() |> 
    mutate(sum_f = sum(Abundance), .by = F) |> 
    mutate(sum_r = sum(Abundance), .by = R) |> 
    mutate(PA_Tissue = ifelse(sum(Abundance) > 0, 1, 0), .by = OTU) |> 
    mutate(fr = log10(sum_f) + log10(sum_r),
           fr_scale = scale(fr)[,1])
  
  phylo_long_colprods
}

correct_phyloseq_tj <- function(ps, ps_long_preds) {
  ## Filter and arrange PS object
  ps_filt <- ps |> 
    ps_filter(SampleType == "Sample", .keep_all_taxa = TRUE) |> 
    ps_arrange(Site, Tree) |>
    tax_sort() |> t()
  
  ## Build matrix!
  thresh_tbl <- ps_long_preds |> 
    ungroup() |> 
    arrange(OTU, Site, Tree) |> 
    select(OTU, Sample, OccThresh) |> 
    pivot_wider(names_from = "Sample", 
                values_from = "OccThresh")
  
  spp <- thresh_tbl$OTU
  
  thresh_mat <- thresh_tbl |> 
    select(-OTU) |> 
    as.matrix()
  rownames(thresh_mat) <- spp
  
  ## Zero things
  tab <- otu_table(ps_filt)
  tab[tab < thresh_mat] <- 0
  
  otu_table(ps_filt) <- tab
  
  ps_filt
}

ps_tagjump_correct <- function(ps) {
  colprods <- estimate_colprods(ps)
  
  mod <-  brms::brm(formula = bf(Abundance ~ 0 + Intercept + fr_scale + (1 + fr_scale | OTU),
                           shape ~ 0 + Intercept + (1 | OTU)),
              family = negbinomial(),
              prior = c(set_prior("normal(0, 3)", class = "b"),
                        set_prior("exponential(2)", class = "sd"),
                        set_prior("lkj(3)", class = "cor")),
              data = dplyr::filter(colprods, 
                                   SampleType == "Unused"),
              chains = 4,
              cores = 4,
              warmup = 1000,
              iter = 3000,
              threads = threading(16),
              backend = "cmdstan",
              silent = 0
  )
  
  ps_long_data <- colprods |>
    filter(SampleType == "Sample") |>
    nest_by(OTU, .keep = TRUE) |>
    mutate(Pred = list(posterior_predict(mod, newdata = data)),
           Quant = list(apply(Pred, 2, \(x) quantile(x, 0.975) |> round())),
           data = list(bind_cols(data, OccThresh = Quant) |> select(-OTU))
    ) |>
    select(-Pred, -Quant) |>
    unnest(cols = c(data))
  
  ps_corrected <- correct_phyloseq_tj(
    ps,
    ps_long_data
  )
  
  return(ps_corrected)
}
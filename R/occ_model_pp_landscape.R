mean_centre_df <- function(input_df, ref_df){
  imap(input_df, \(x, idx) (x - mean(pull(ref_df, any_of(idx)))) / sd(pull(ref_df, any_of(idx)))) |> 
    as_tibble()
}

landscape_predictions <- function(draws, 
                                  ref_df, 
                                  landscape_df, 
                                  sample_df, 
                                  taxa = NULL,
                                  chunk_size = NULL,
                                  .fun_pp = posterior_predict_phi,
                                  phi = NULL,
                                  theta = NULL){
  if(!all(names(landscape_df) %in% names(ref_df))) {
    stop("Error: missing variable in input data frame!")
  }
  
  draws <- posterior::as_draws_rvars(draws)
  n_taxa <- dim(draws$b_env_species)[3]
  
  if(is.null(taxa)) {
    taxa <- 1:n_taxa
  }
  
  if(is.null(chunk_size)) {
    taxa_chunks <- taxa
  } else {
    z <- seq_along(taxa)
    taxa_chunks <- split(taxa, ceiling(z/chunk_size))
    names(taxa_chunks) <- NULL
  }
  
  landscape_scaled <- select(landscape_df, -Easting, -Northing) |> 
    mean_centre_df(ref_df)

  env_modmat <- model.matrix(~ ., data = landscape_scaled)
  sample_modmat <- model.matrix(~ 0 + ., data = sample_df)

  map(taxa_chunks, \(x) .fun_pp(x, draws, env_modmat, sample_modmat, phi, theta), .progress = TRUE)
}

bind_pred_landscape <- function(landscape_df, prediction) {
  select(landscape_df, Easting, Northing) |> 
    mutate(pred = as.vector(prediction))
}

plot_raster <- function(predictions, landscape_df, site_df, ps, asv_tbl, asv_index, .predfun = mean) {
  asv <- asv_tbl |> 
    filter(ASV_numeric == asv_index) |>
    mutate(across(everything(), \(x) ifelse(x == "", "unknown", x)))
  
  ps2 <- speedyseq::psmelt(ps) |> 
    ungroup() |> 
    select(OTU, Site, Abundance) |> 
    summarise(Pres = any(Abundance > 0) * 1, 
              Prop = sum(Abundance > 0) / n(),
              .by = c(OTU, Site)) |> 
    arrange(OTU, Site) |> 
    mutate(ASV_numeric = as.numeric(as.factor(OTU))) |> 
    filter(ASV_numeric == asv_index)
    
  site_df <- site_df |>
    sf::st_as_sf(coords = c("Easting", "Northing"), crs = 27700) |>
    left_join(ps2)

  mutate(landscape_df, pred = .predfun(predictions[,asv_index])) |> 
    ggplot() +
    geom_tile(aes(x = Easting, y = Northing, fill = pred)) +
    geom_sf(data = site_df, aes(size = Prop, shape = as.factor(Pres)), colour = "white") +
    coord_sf(crs = 27700) +
    scale_fill_viridis_c(na.value = 0, limits = c(0, 1)) + 
    labs(title = glue::glue_data(asv, "ASV ID: {ASV_ID}
    Phylum: {Phylum}
    Class: {Class}
    Order: {Order}
    Family: {Family}
    Binomial: {Genus} {Species}"),
         fill = "probability taxon present",
         x = "", y = "") +
    theme_light() +
    scale_shape_manual(values = c(1, 16)) +
    scale_size_continuous(limits = c(0, 1)) +
    guides(shape = "none",
           fill = guide_colourbar(title.position="top")) +
    theme(legend.position = "right",
          plot.title = element_text(size = 8))
}

write_maps <- function(map, taxon, path, asv_table) {
  if(!dir.exists(path)) dir.create(path, recursive = TRUE)
  taxid <- asv_table |> filter(ASV_numeric == taxa) |> pull(ASV_ID)
  ggsave(plot = map, filename = paste0(dir, "/", taxid, ".png"),
         width = 6, height = 7)
}

# spat <- mutate(geom, pred = preds$V1) |>
#   mutate(Easting = sf::st_coordinates(geom)[,1],
#          Northing = sf::st_coordinates(geom)[,2]) |>
#   sf::st_drop_geometry() |>
#   select(Easting, Northing, pred) |>
#   mutate(Easting = round(Easting),
#          Northing = round(Northing)) |> 
#   tidyterra::as_spatraster(crs = "EPSG:27700", digits = 3) |> 
#   terra::wrap()
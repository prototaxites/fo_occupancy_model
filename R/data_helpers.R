## Filter a community count matrix to remove rows with zero total reads
ps_filter_rowsum <- function(ps, min_total = 0) {
  ps_filtered <- ps |> 
    ps_mutate(temp_N = ifelse(
      taxa_are_rows(ps),
      colSums(otu_table(ps)),
      rowSums(otu_table(ps))
    )) |> 
    phyloseq::ps_filter(temp_N > 0) |> 
    select(-temp_N)
  
  return(ps_filtered)
}

## Filter a community count matrix to remove ASVs
## that do not make up some minimum proportion of reads in at
## least one sample.
ps_filter_minfrac <- function(ps, min_proportion = 0.05) {
  keep_asvs <- ps |> 
    transform_sample_counts(\(x) x / sum(x)) |> 
    phyloseq::filter_taxa(\(x) any(x > min_proportion))
  
  ps_filtered <- ps |>
    prune_taxa(keep_asvs, x = _)
  
  return(ps_filtered)
}



read_landscape <- function(landscape_file, vars) {
  df <- sf::read_sf(landscape_file, as_tibble = TRUE) |> 
    rename(groundfrost = grndfrs, 
           rainfall = rainfll, 
           snowLying = snwLyng) |> 
    select(all_of(vars))
  
  idx <- complete.cases(sf::st_drop_geometry(df))
  
  filter(df, idx) |> 
    mutate(Easting = sf::st_coordinates(geometry)[,1],
           Northing = sf::st_coordinates(geometry)[,2]) |> 
    sf::st_drop_geometry()
}

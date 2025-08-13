plot_disease_intervals <- function(disease_intervals) {
  disease_intervals |> 
    mutate(n = 1:n()) |> 
    mutate(tax = glue::glue("{Kingdom}: {Phylum}")) |> 
    ggplot(aes(y = tax, colour = Tissue)) +
    # geom_rect(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#FF7256", alpha = 0.005) +
    # geom_rect(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#C6E2FF", alpha = 0.01) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 1) +
    stat_pointinterval(aes(xdist = diff_lfc), .width = 0.95, linewidth = 0.2, position = position_jitter(width = 0.05)) +
    scale_colour_manual(values = c(Leaf = "olivedrab3", Stem = "bisque2", Soil = "salmon4")) +
    labs(y = NULL, 
         x = "log2 fold change in read abundance") +
    theme_bw() +
    theme(#axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

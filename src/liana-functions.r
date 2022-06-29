#' Customization of liana_dotplot() function
#'
#' @details We modified color palette of axis.text.y to
#' scales::hue_pal(length(unique(liana_mod$target))) and element text sizes.
#' For general information on liana_doplot()
#' see https://rdrr.io/github/saezlab/ligrec_decoupler/src/R/liana_plot.R
my_liana_dotplot <- function(liana_agg, source_groups, target_groups, specificity = "natmi.edge_specificity",
                             magnitude = "sca.LRscore", show_complex = TRUE) {
  if (show_complex) {
    entities <- c("ligand.complex", "receptor.complex")
  } else {
    entities <- c("ligand", "receptor")
  }
  liana_mod <- liana_agg %>%
    filter(source %in% source_groups) %>%
    filter(target %in% target_groups) %>%
    rename(magnitude = !!magnitude) %>%
    rename(specificity = !!specificity) %>%
    unite(entities,
      col = "interaction", sep = " -> "
    ) %>%
    unite(c(
      "source",
      "target"
    ), col = "source_target", remove = FALSE)
  cbPalette <- scales::hue_pal()(20)
  suppressWarnings(ggplot(liana_mod, aes(
    x = interaction, y = target,
    colour = magnitude, size = specificity, group = target
  )) +
    geom_point() +
    scale_color_gradientn(colours = viridis::viridis(20)) +
    scale_size_continuous(range = c(5, 9)) +
    facet_grid(source ~
      ., space = "free", scales = "free", switch = "y") +
    theme_bw(base_size = 20) +
    theme(
      legend.text = element_text(size = 12), axis.text.y = element_text(
        colour = scales::hue_pal()(length(unique(liana_mod$target))),
        face = "bold", size = 12
      ), axis.text.x = element_text(
        size = 12,
        angle = 90, vjust = 0.5
      ), legend.title = element_text(size = 12),
      panel.spacing = unit(0.1, "lines"), strip.background = element_rect(fill = NA),
      strip.text = element_text(size = 5, colour = "gray6")
    ) +
    scale_y_discrete(position = "right") +
    labs(
      x = "Interactions (Ligand -> Receptor)",
      colour = "Expression\nMagnitude", size = "Interaction\nSpecificity",
      y = NULL
    ))
}

#' Customization of liana_dotplot() function
#'
#' @details We modified color palette of axis.text.y to "black"
#' and element text sizes. For general information on liana_doplot()
#' see https://rdrr.io/github/saezlab/ligrec_decoupler/src/R/liana_plot.R
my_liana_dotplot2 <- function(liana_agg, source_groups, target_groups, specificity = "natmi.edge_specificity",
                              magnitude = "sca.LRscore", show_complex = TRUE) {
  if (show_complex) {
    entities <- c("ligand.complex", "receptor.complex")
  } else {
    entities <- c("ligand", "receptor")
  }
  liana_mod <- liana_agg %>%
    filter(source %in% source_groups) %>%
    filter(target %in% target_groups) %>%
    rename(magnitude = !!magnitude) %>%
    rename(specificity = !!specificity) %>%
    unite(entities,
      col = "interaction", sep = " -> "
    ) %>%
    unite(c(
      "source",
      "target"
    ), col = "source_target", remove = FALSE)
  cbPalette <- scales::hue_pal()(20)
  suppressWarnings(ggplot(liana_mod, aes(
    x = interaction, y = target,
    colour = magnitude, size = specificity, group = target
  )) +
    geom_point() +
    scale_color_gradientn(colours = viridis::viridis(20)) +
    scale_size_continuous(range = c(5, 9)) +
    facet_grid(source ~
      ., space = "free", scales = "free", switch = "y") +
    theme_bw(base_size = 20) +
    theme(
      legend.text = element_text(size = 16), axis.text.y = element_text(
        colour = "black",
        face = "bold", size = 23
      ), axis.text.x = element_text(
        size = 18,
        angle = 90, vjust = 0.5
      ), legend.title = element_text(size = 18),
      panel.spacing = unit(0.1, "lines"), strip.background = element_rect(fill = NA),
      strip.text = element_text(size = 16, colour = "gray6")
    ) +
    scale_y_discrete(position = "right") +
    labs(
      x = "Interactions (Ligand -> Receptor)",
      colour = "Expression\nMagnitude", size = "Interaction\nSpecificity",
      y = NULL
    ))
}

plot_lymphs_vs_mels <- function(liana_obj, seurat_obj, dataset) {
  lymphocytes <- c("Tc", "Th", "Treg", "NK1", "NK2", "NK3", "Plasma cell", "ILC")
  tcells <- c("Tc", "Th", "Treg")
  melanocytes <- levels(Idents(seurat_obj))[grep("Melanocyte", levels(Idents(seurat_obj)))]

  filename <- here(figures_path, paste0("06_", dataset, "-rl-analysis-lymphs-vs-mels.pdf"))
  pdf(filename, width = 28, height = 14)

  liana_subset <- liana_obj %>%
    # filter(source == melanocytes ) %>%
    filter(aggregate_rank < 0.05)
  if (nrow(liana_subset) > 0) {
    my_liana_dotplot2(liana_subset,
      source_groups = melanocytes,
      target_groups = lymphocytes
    ) %>%
      print() %>%
      tryCatch()
  }

  liana_subset <- liana_obj %>%
    # filter(source == lymphocytes ) %>%
    filter(aggregate_rank < 0.05)
  if (nrow(liana_subset) > 0) {
    my_liana_dotplot2(liana_subset,
      source_groups = lymphocytes,
      target_groups = melanocytes
    ) %>%
      print() %>%
      tryCatch()
  }


  liana_subset <- liana_obj %>%
    # filter(source == melanocytes ) %>%
    filter(aggregate_rank < 0.05)
  if (nrow(liana_subset) > 0) {
    my_liana_dotplot2(liana_subset,
      source_groups = melanocytes,
      target_groups = tcells
    ) %>%
      print() %>%
      tryCatch()
  }

  liana_subset <- liana_obj %>%
    # filter(source == tcells ) %>%
    filter(aggregate_rank < 0.05)
  if (nrow(liana_subset) > 0) {
    my_liana_dotplot2(liana_subset,
      source_groups = tcells,
      target_groups = melanocytes
    ) %>%
      print() %>%
      tryCatch()
  }

  dev.off()
}

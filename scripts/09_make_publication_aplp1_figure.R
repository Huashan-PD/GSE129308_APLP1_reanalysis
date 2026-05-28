#!/usr/bin/env Rscript

source("scripts/lib_project.R")

run_main("09_make_publication_aplp1_figure", {
  load_package("data.table", required = TRUE, log = "09_make_publication_aplp1_figure.log")
  load_package("ggplot2", required = TRUE, log = "09_make_publication_aplp1_figure.log")
  load_package("patchwork", required = TRUE, log = "09_make_publication_aplp1_figure.log")
  load_package("scales", required = TRUE, log = "09_make_publication_aplp1_figure.log")

  summary <- data.table::fread(project_path("outputs/tables/aplp1_summary_by_condition_subtype.csv"))
  fraction <- data.table::fread(project_path("outputs/tables/aplp1_fraction_by_donor_condition_subtype.csv"))
  stats_all <- data.table::fread(project_path("outputs/tables/aplp1_stats_all.csv"))

  subtype_order <- c(
    "Ex1", "Ex2", "Ex3", "Ex4/5", "Ex6", "Ex7", "Ex8", "Ex10", "Ex11", "Ex12", "Ex13",
    "In1", "In2", "In3/4", "In5/6", "In7",
    "Excitatory_unresolved", "Inhibitory_unresolved", "Unassigned"
  )
  observed_subtypes <- intersect(subtype_order, unique(summary$neuronal_subtype))
  observed_subtypes <- c(observed_subtypes, setdiff(sort(unique(summary$neuronal_subtype)), observed_subtypes))

  summary[, condition := factor(condition, levels = condition_levels)]
  summary[, neuronal_subtype := factor(neuronal_subtype, levels = rev(observed_subtypes))]
  fraction[, condition := factor(condition, levels = condition_levels)]
  fraction[, neuronal_subtype := factor(neuronal_subtype, levels = rev(observed_subtypes))]

  source_data <- summary[
    ,
    .(
      neuronal_subtype,
      condition,
      n_donors,
      total_cells,
      total_APLP1_positive,
      pooled_APLP1_positive_fraction,
      median_donor_APLP1_positive_fraction,
      mean_donor_normalized_expression
    )
  ]
  write_table(source_data, project_path("outputs/tables/Fig_APLP1_publication_source_data.csv"))

  plot_theme <- ggplot2::theme_classic(base_size = 7, base_family = "Helvetica") +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.25),
      legend.title = ggplot2::element_text(size = 6.5),
      legend.text = ggplot2::element_text(size = 6.5),
      plot.title = ggplot2::element_text(size = 7.2, hjust = 0.5, face = "plain", lineheight = 0.95),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )

  heat_scale <- ggplot2::scale_fill_gradientn(
    colors = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    limits = c(0.25, 0.80),
    oob = scales::squish,
    breaks = c(0.30, 0.40, 0.50, 0.60, 0.70),
    name = "Median donor\nAPLP1-positive\nfraction"
  )

  example_dt <- summary[condition %in% c("AD_AT8pos", "nonAD")]
  example_dt[, condition_plot := factor(
    as.character(condition),
    levels = c("AD_AT8pos", "nonAD"),
    labels = c("AD-AT8+", "non-AD\ncontrol")
  )]

  p_example <- ggplot2::ggplot(
    example_dt,
    ggplot2::aes(x = condition_plot, y = neuronal_subtype, fill = median_donor_APLP1_positive_fraction)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.45, width = 0.95, height = 0.95) +
    heat_scale +
    ggplot2::coord_fixed(ratio = 0.55) +
    ggplot2::labs(
      title = "Fraction of AD-AT8+ and non-AD control\nneuronal somas expressing APLP1",
      x = NULL,
      y = NULL
    ) +
    plot_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 7, face = "bold"),
      legend.key.height = grid::unit(0.30, "cm"),
      legend.key.width = grid::unit(0.20, "cm")
    )

  ggplot2::ggsave(
    project_path("outputs/figures/Fig_APLP1_fraction_heatmap_example_style.pdf"),
    p_example,
    width = 3.25,
    height = 2.65,
    useDingbats = FALSE
  )
  ggplot2::ggsave(
    project_path("outputs/figures/Fig_APLP1_fraction_heatmap_example_style.png"),
    p_example,
    width = 3.25,
    height = 2.65,
    dpi = 600
  )

  primary_dt <- summary[condition %in% c("AD_AT8neg", "AD_AT8pos")]
  primary_dt[, condition_plot := factor(
    as.character(condition),
    levels = c("AD_AT8neg", "AD_AT8pos"),
    labels = c("AD-\nAT8-", "AD-\nAT8+")
  )]
  p_primary_heat <- ggplot2::ggplot(
    primary_dt,
    ggplot2::aes(x = condition_plot, y = neuronal_subtype, fill = median_donor_APLP1_positive_fraction)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.45, width = 0.95, height = 0.95) +
    heat_scale +
    ggplot2::coord_fixed(ratio = 0.55) +
    ggplot2::labs(title = "Primary fraction", x = NULL, y = NULL) +
    plot_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 7, face = "bold"),
      legend.position = "none"
    )

  secondary_dt <- example_dt
  secondary_dt[, condition_plot := factor(
    as.character(condition),
    levels = c("AD_AT8pos", "nonAD"),
    labels = c("AD-\nAT8+", "non-AD\ncontrol")
  )]
  p_secondary_heat <- ggplot2::ggplot(
    secondary_dt,
    ggplot2::aes(x = condition_plot, y = neuronal_subtype, fill = median_donor_APLP1_positive_fraction)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.45, width = 0.95, height = 0.95) +
    heat_scale +
    ggplot2::coord_fixed(ratio = 0.55) +
    ggplot2::labs(title = "AD-AT8+ vs non-AD", x = NULL, y = NULL) +
    plot_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7, face = "bold"),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "right"
    )

  pb <- stats_all[
    comparison == "AD_AT8pos_vs_AD_AT8neg" &
      metric == "APLP1_pseudobulk_expression"
  ]
  pb[, neuronal_subtype := factor(neuronal_subtype, levels = rev(observed_subtypes))]
  pb[, sig := data.table::fcase(
    status != "ok", "insufficient",
    !is.na(FDR) & FDR < 0.05, "FDR < 0.05",
    default = "n.s."
  )]
  p_pb <- ggplot2::ggplot(pb, ggplot2::aes(x = log2FC, y = neuronal_subtype)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey55") +
    ggplot2::geom_point(ggplot2::aes(color = sig), size = 1.8, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c("FDR < 0.05" = "#b2182b", "n.s." = "#4d4d4d", "insufficient" = "#bdbdbd"),
      breaks = c("FDR < 0.05", "n.s.", "insufficient"),
      name = NULL
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "Primary pseudobulk expression",
      x = "log2FC (AD-AT8+ / AD-AT8-)",
      y = NULL
    ) +
    plot_theme +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7, face = "bold"),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )

  p_multi <- (p_primary_heat | p_secondary_heat | p_pb) +
    patchwork::plot_layout(widths = c(0.95, 1.15, 1.75), guides = "keep") +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 9))

  ggplot2::ggsave(
    project_path("outputs/figures/Fig_APLP1_publication_multipanel.pdf"),
    p_multi,
    width = 8.4,
    height = 3.25,
    useDingbats = FALSE
  )
  ggplot2::ggsave(
    project_path("outputs/figures/Fig_APLP1_publication_multipanel.png"),
    p_multi,
    width = 8.4,
    height = 3.25,
    dpi = 600
  )

  legend_text <- c(
    "Figure. APLP1 detection and expression in GSE129308 neuronal somas.",
    "",
    "A, Heatmap showing the median donor-level APLP1-positive fraction, defined as the fraction of neuronal somas with raw APLP1 counts > 0, in paired AD-AT8- neighboring NFT-free and AD-AT8+ NFT-bearing neuronal somas. B, Example-style heatmap comparing AD-AT8+ neuronal somas with non-AD control neuronal somas. non-AD controls are not labeled as AT8-negative because AT8 status was defined by FACS immunolabeling in the AD samples, not inferred from RNA-seq. C, Donor-level pseudobulk APLP1 expression in the primary paired comparison, computed by aggregating counts to donor x condition x neuronal subtype and testing AD-AT8+ versus AD-AT8- with limma-voom/edgeR using design ~ donor + condition. Points show log2 fold change. Fraction statistics were tested at the donor level with paired Wilcoxon tests for AD-AT8+ versus AD-AT8- and Wilcoxon rank-sum tests for AD groups versus non-AD controls; p values were adjusted by Benjamini-Hochberg FDR across neuronal subtypes. Subtypes not meeting the prespecified donor/cell-count threshold are labeled insufficient. The figure distinguishes APLP1-positive fraction from APLP1 expression level; fraction differences alone should be described as changes in the APLP1-positive fraction, not expression upregulation."
  )
  writeLines(legend_text, project_path("outputs/report/Fig_APLP1_publication_legend.txt"))

  log_message("Publication APLP1 figures written to outputs/figures.", log = "09_make_publication_aplp1_figure.log")
})

#!/usr/bin/env Rscript

source("scripts/lib_project.R")

run_main("11_make_supplementary_s2_aplp1", {
  load_package("data.table", required = TRUE, log = "11_make_supplementary_s2_aplp1.log")
  load_package("ggplot2", required = TRUE, log = "11_make_supplementary_s2_aplp1.log")
  load_package("patchwork", required = TRUE, log = "11_make_supplementary_s2_aplp1.log")
  load_package("scales", required = TRUE, log = "11_make_supplementary_s2_aplp1.log")

  author_source_path <- project_path("outputs/tables/author_tableS6_APLP1_fraction_author_subtypes_3groups.csv")
  donor_fraction_path <- project_path("outputs/tables/aplp1_fraction_by_donor_condition_subtype.csv")
  stats_path <- project_path("outputs/tables/aplp1_stats_all.csv")
  for (p in c(author_source_path, donor_fraction_path, stats_path)) {
    if (!file.exists(p)) stop("Required input is missing: ", p, call. = FALSE)
  }

  row_order <- c("Ex1", "Ex2", "Ex3", "Ex4/5", "Ex6", "Ex7", "Ex8", "Ex10", "Ex13", "In2")
  row_order_author <- gsub("/", "&", row_order)
  title_detectable <- "Fraction of neuronal somas with detectable APLP1 transcripts\nin the Otero-Garcia BA9 single-soma dataset"

  format_num <- function(x, digits = 3) {
    ifelse(is.na(x), "n/a", formatC(x, digits = digits, format = "f"))
  }
  format_p <- function(x) {
    ifelse(is.na(x), "n/a", ifelse(x < 0.001, formatC(x, format = "e", digits = 2), formatC(x, format = "f", digits = 3)))
  }

  author_dt <- data.table::fread(author_source_path)
  author_dt[, neuronal_subtype := gsub("&", "/", neuronal_subtype, fixed = TRUE)]
  author_dt[, condition := factor(condition, levels = c("non-AD", "AD-AT8-", "AD-AT8+"))]
  author_dt[, neuronal_subtype := factor(neuronal_subtype, levels = rev(row_order))]
  author_dt[, missing_reason := data.table::fcase(
    is.na(APLP1_positive_fraction) & !is.na(author_n_cells), "not available",
    is.na(APLP1_positive_fraction) & is.na(author_n_cells), "insufficient / not available",
    default = NA_character_
  )]
  write_table(author_dt, project_path("outputs/tables/SuppFig_S2A_author_heatmap_source.csv"))

  p_s2a <- ggplot2::ggplot(
    author_dt,
    ggplot2::aes(x = condition, y = neuronal_subtype, fill = APLP1_positive_fraction)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35, width = 0.96, height = 0.96) +
    ggplot2::geom_text(
      data = author_dt[is.na(APLP1_positive_fraction)],
      ggplot2::aes(label = "n/a"),
      color = "#555555",
      size = 2.0,
      inherit.aes = TRUE
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c("#f7f4ff", "#c5baff", "#7567ff", "#2b00d4"),
      limits = c(0.28, 0.70),
      breaks = c(0.30, 0.40, 0.50, 0.60, 0.70),
      labels = scales::number_format(accuracy = 0.1),
      oob = scales::squish,
      na.value = "#eeeeee",
      name = "APLP1-positive\nfraction"
    ) +
    ggplot2::scale_x_discrete(labels = c("non-AD", "AD-AT8-", "AD-AT8+")) +
    ggplot2::coord_fixed(ratio = 0.30) +
    ggplot2::labs(title = "A  Author subtype heatmap", x = NULL, y = NULL) +
    ggplot2::theme_classic(base_size = 7, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 8, face = "bold", color = "black"),
      axis.text.x = ggplot2::element_text(size = 7, color = "black", margin = ggplot2::margin(t = 2)),
      axis.text.y = ggplot2::element_text(size = 7, face = "bold", color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.25),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.25),
      legend.key.height = grid::unit(0.36, "cm"),
      legend.key.width = grid::unit(0.22, "cm"),
      legend.text = ggplot2::element_text(size = 6.3, color = "black"),
      legend.title = ggplot2::element_text(size = 6.5, color = "black"),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  key_subtypes <- c("Ex2", "Ex7", "Ex10")
  donor_fraction <- data.table::fread(donor_fraction_path)
  paired_dt <- donor_fraction[
    neuronal_subtype %in% key_subtypes &
      condition %in% c("AD_AT8neg", "AD_AT8pos") &
      n_cells >= 30
  ]
  paired_wide <- data.table::dcast(
    paired_dt,
    donor + neuronal_subtype ~ condition,
    value.var = c("APLP1_positive_fraction", "n_cells")
  )
  paired_wide <- paired_wide[
    !is.na(APLP1_positive_fraction_AD_AT8neg) &
      !is.na(APLP1_positive_fraction_AD_AT8pos)
  ]
  paired_source <- data.table::melt(
    paired_wide,
    id.vars = c("donor", "neuronal_subtype"),
    measure.vars = c("APLP1_positive_fraction_AD_AT8neg", "APLP1_positive_fraction_AD_AT8pos"),
    variable.name = "condition_source",
    value.name = "APLP1_positive_fraction"
  )
  paired_source[, condition := sub("APLP1_positive_fraction_", "", condition_source)]
  paired_source[, condition := factor(condition, levels = c("AD_AT8neg", "AD_AT8pos"), labels = c("AD-AT8-", "AD-AT8+"))]
  paired_source[, neuronal_subtype := factor(neuronal_subtype, levels = key_subtypes)]
  write_table(paired_source, project_path("outputs/tables/SuppFig_S2B_paired_donor_source.csv"))

  p_s2b <- ggplot2::ggplot(
    paired_source,
    ggplot2::aes(x = condition, y = APLP1_positive_fraction, group = donor)
  ) +
    ggplot2::geom_line(color = "#777777", linewidth = 0.35, alpha = 0.85) +
    ggplot2::geom_point(ggplot2::aes(fill = condition), shape = 21, size = 1.7, color = "black", stroke = 0.25) +
    ggplot2::facet_wrap(~ neuronal_subtype, nrow = 1) +
    ggplot2::scale_fill_manual(values = c("AD-AT8-" = "#a6bddb", "AD-AT8+" = "#2b00d4"), guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0.25, 0.90)) +
    ggplot2::labs(
      title = "B  Paired donor-level fractions",
      x = NULL,
      y = "Detectable APLP1 fraction"
    ) +
    ggplot2::theme_classic(base_size = 7, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 8, face = "bold", color = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 7, face = "bold", color = "black"),
      axis.text.x = ggplot2::element_text(size = 7, color = "black"),
      axis.text.y = ggplot2::element_text(size = 6.5, color = "black"),
      axis.title.y = ggplot2::element_text(size = 7, color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.25),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.25),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  stats <- data.table::fread(stats_path)
  frac_stats <- stats[
    comparison == "AD_AT8pos_vs_AD_AT8neg" &
      metric == "APLP1_positive_fraction"
  ]
  pb_stats <- stats[
    comparison == "AD_AT8pos_vs_AD_AT8neg" &
      metric == "APLP1_pseudobulk_expression",
    .(neuronal_subtype, pseudobulk_log2FC = log2FC, pseudobulk_FDR = FDR, pseudobulk_status = status)
  ]
  n_cells <- donor_fraction[
    condition %in% c("AD_AT8neg", "AD_AT8pos"),
    .(n_cells_total = sum(n_cells), n_donors = uniqueN(donor[n_cells >= 30])),
    by = .(neuronal_subtype, condition)
  ]
  n_cells_wide <- data.table::dcast(
    n_cells,
    neuronal_subtype ~ condition,
    value.var = c("n_cells_total", "n_donors")
  )
  stat_summary <- merge(frac_stats, pb_stats, by = "neuronal_subtype", all = TRUE)
  stat_summary <- merge(stat_summary, n_cells_wide, by = "neuronal_subtype", all.x = TRUE)
  stat_summary <- stat_summary[
    neuronal_subtype %in% c("Ex2", "Ex7", "Ex10", "Ex1", "In5/6", "In1", "Ex6", "Ex11", "In3/4")
  ]
  stat_summary[, display_subtype := neuronal_subtype]
  stat_summary[, n_donors_display := ifelse(is.na(n_paired_donors), "insuff.", as.character(n_paired_donors))]
  stat_summary[, n_cells_display := ifelse(
    is.na(n_cells_total_AD_AT8neg) | is.na(n_cells_total_AD_AT8pos),
    "n/a",
    paste0(scales::comma(n_cells_total_AD_AT8neg), "/", scales::comma(n_cells_total_AD_AT8pos))
  )]
  stat_summary[, delta_fraction_display := format_num(effect_size, 3)]
  stat_summary[, log2FC_display := format_num(pseudobulk_log2FC, 3)]
  stat_summary[, FDR_display := format_p(FDR)]
  stat_summary[, pb_FDR_display := format_p(pseudobulk_FDR)]
  stat_summary[, status_display := ifelse(status == "ok" & pseudobulk_status == "ok", "ok", "limited")]
  stat_summary[, sort_i := data.table::fcase(
    neuronal_subtype == "Ex2", 1L,
    neuronal_subtype == "Ex7", 2L,
    neuronal_subtype == "Ex10", 3L,
    neuronal_subtype == "Ex1", 4L,
    neuronal_subtype == "In5/6", 5L,
    neuronal_subtype == "In1", 6L,
    neuronal_subtype == "Ex6", 7L,
    neuronal_subtype == "Ex11", 8L,
    neuronal_subtype == "In3/4", 9L,
    default = 99L
  )]
  data.table::setorder(stat_summary, sort_i)
  write_table(stat_summary, project_path("outputs/tables/SuppFig_S2C_primary_stats_summary.csv"))

  table_dt <- stat_summary[
    ,
    .(
      row_i = .I,
      subtype = display_subtype,
      donors = n_donors_display,
      cells = n_cells_display,
      delta = delta_fraction_display,
      log2FC = log2FC_display,
      fraction_FDR = FDR_display,
      expression_FDR = pb_FDR_display
    )
  ]
  table_long <- data.table::melt(
    table_dt,
    id.vars = "row_i",
    variable.name = "column",
    value.name = "value"
  )
  column_levels <- c("subtype", "donors", "cells", "delta", "log2FC", "fraction_FDR", "expression_FDR")
  column_labels <- c(
    subtype = "Subtype",
    donors = "Paired\ndonors",
    cells = "Cells\nAT8-/AT8+",
    delta = "Delta\nfraction",
    log2FC = "Pseudo-\nbulk log2FC",
    fraction_FDR = "Fraction\nFDR",
    expression_FDR = "Expr.\nFDR"
  )
  table_long[, column := factor(column, levels = column_levels, labels = column_labels[column_levels])]
  table_long[, y := max(row_i) - row_i + 1]
  header_dt <- data.table::data.table(
    column = factor(column_labels[column_levels], levels = column_labels[column_levels]),
    y = max(table_long$y) + 1,
    value = unname(column_labels[column_levels])
  )

  p_s2c <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = table_long,
      ggplot2::aes(x = column, y = y),
      fill = "#f7f7f7",
      color = "white",
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      data = table_long,
      ggplot2::aes(x = column, y = y, label = value),
      size = 2.1,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = header_dt,
      ggplot2::aes(x = column, y = y, label = value),
      size = 2.0,
      fontface = "bold",
      color = "black",
      lineheight = 0.9
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
    ggplot2::labs(title = "C  Primary paired statistics", x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 7, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 8, face = "bold", color = "black", hjust = 0),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  save_panel <- function(plot, name, width, height) {
    ggplot2::ggsave(project_path("outputs/figures", paste0(name, ".pdf")), plot, width = width, height = height, useDingbats = FALSE)
    ggplot2::ggsave(project_path("outputs/figures", paste0(name, ".png")), plot, width = width, height = height, dpi = 600)
  }

  save_panel(p_s2a, "SuppFig_S2A_author_subtype_three_group_heatmap", 3.65, 2.45)
  save_panel(p_s2b, "SuppFig_S2B_paired_donor_key_subtypes", 4.3, 2.25)
  save_panel(p_s2c, "SuppFig_S2C_primary_stats_summary", 6.4, 2.35)

  top_row <- p_s2a | p_s2b
  supp_s2 <- top_row / p_s2c +
    patchwork::plot_layout(heights = c(1.0, 0.92)) +
    patchwork::plot_annotation(
      title = title_detectable,
      caption = paste(
        "S2A uses author Table S6 subtype-level pct values; grey/n/a indicates not available in the APLP1 row.",
        "S2B-S2C use donor-level H5 reanalysis, treating donor rather than cell as the replicate.",
        "Delta fraction and pseudobulk log2FC are AD-AT8+ minus/over AD-AT8-."
      )
    ) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 9, color = "black", lineheight = 0.95),
      plot.caption = ggplot2::element_text(size = 5.8, color = "#333333", hjust = 0)
    )

  save_panel(supp_s2, "SuppFig_S2_APLP1_detectable_fraction", 8.4, 5.2)

  legend_text <- c(
    "Supplementary Fig. S2. Detectable APLP1 transcripts in neuronal somas from the Otero-Garcia BA9 single-soma dataset.",
    "",
    "A, Heatmap of author Table S6 APLP1-positive fractions across non-AD, AD-AT8-, and AD-AT8+ groups using the author neuronal subtype order. APLP1-positive fraction denotes raw APLP1 UMI counts > 0. Grey/n/a indicates that the value was not available in the author Table S6 APLP1 row. B, Paired donor-level APLP1-positive fractions for selected reanalyzed neuronal subtypes, comparing AD-AT8- neighboring NFT-free neuronal somas with matched AD-AT8+ NFT-bearing neuronal somas. Each line is one AD donor. C, Primary paired comparison summary from donor-level reanalysis. Fraction tests use donor-level paired comparisons, and pseudobulk expression aggregates counts to donor x condition x subtype before testing. non-AD controls are non-AD control neuronal somas and are not treated as an AT8-negative group."
  )
  writeLines(legend_text, project_path("outputs/report/SuppFig_S2_APLP1_detectable_fraction_legend.txt"))

  log_message("Supplementary Fig. S2 written to outputs/figures.", log = "11_make_supplementary_s2_aplp1.log")
})

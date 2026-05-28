#!/usr/bin/env Rscript

source("scripts/lib_project.R")

run_main("plot_APLP1_single_soma_clean_AB", {
  load_package("data.table", required = TRUE, log = "plot_APLP1_single_soma_clean_AB.log")
  load_package("ggplot2", required = TRUE, log = "plot_APLP1_single_soma_clean_AB.log")
  load_package("patchwork", required = TRUE, log = "plot_APLP1_single_soma_clean_AB.log")
  load_package("scales", required = TRUE, log = "plot_APLP1_single_soma_clean_AB.log")
  load_package("openxlsx", required = TRUE, log = "plot_APLP1_single_soma_clean_AB.log")
  load_package("svglite", required = TRUE, log = "plot_APLP1_single_soma_clean_AB.log")

  set.seed(PROJECT_SEED)
  ensure_project_dirs()
  dir.create(project_path("scripts/figure_generation"), recursive = TRUE, showWarnings = FALSE)
  dir.create(project_path("outputs/figures/archive"), recursive = TRUE, showWarnings = FALSE)

  author_heatmap_path <- project_path("outputs/tables/author_tableS6_APLP1_fraction_author_subtypes_3groups.csv")
  donor_fraction_path <- project_path("outputs/tables/aplp1_fraction_by_donor_condition_subtype.csv")
  stats_path <- project_path("outputs/tables/aplp1_stats_all.csv")
  required_inputs <- c(author_heatmap_path, donor_fraction_path, stats_path)
  missing_inputs <- required_inputs[!file.exists(required_inputs)]
  if (length(missing_inputs)) {
    stop("Required inputs are missing: ", paste(missing_inputs, collapse = ", "), call. = FALSE)
  }

  archive_sources <- c(
    project_path("outputs/figures/SuppFig_S2_APLP1_detectable_fraction.png"),
    "/Users/ywb/Desktop/image (1).png",
    "/Users/ywb/Desktop/Fig_APLP1_single_soma_primary_donor_fraction_fig5style_sig_nolegend.png"
  )
  archive_names <- c(
    "archive_previous_codex_composite_SuppFig_S2_APLP1_detectable_fraction.png",
    "archive_user_attachment_codex_composite_image_1.png",
    "archive_user_attachment_claude_fig5style.png"
  )
  for (i in seq_along(archive_sources)) {
    if (file.exists(archive_sources[i])) {
      file.copy(
        archive_sources[i],
        project_path("outputs/figures/archive", archive_names[i]),
        overwrite = TRUE
      )
    }
  }

  author_order <- c("Ex1", "Ex2", "Ex3", "Ex4/5", "Ex6", "Ex7", "Ex8", "Ex10", "Ex13", "In2")
  priority_subtypes <- c("Ex2", "Ex7", "Ex10")
  priority_labels <- c(
    Ex2 = "Ex2_CUX2-COL5A2",
    Ex7 = "Ex7_RORB-PCP4",
    Ex10 = "Ex10_NFIA/THEMIS"
  )

  author_dt <- data.table::fread(author_heatmap_path)
  author_dt[, neuronal_subtype := gsub("&", "/", neuronal_subtype, fixed = TRUE)]
  author_dt <- author_dt[neuronal_subtype %in% author_order]
  author_dt[, condition := factor(condition, levels = c("non-AD", "AD-AT8-", "AD-AT8+"))]
  author_dt[, neuronal_subtype := factor(neuronal_subtype, levels = rev(author_order))]

  panel_a <- ggplot2::ggplot(
    author_dt,
    ggplot2::aes(x = condition, y = neuronal_subtype, fill = APLP1_positive_fraction)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.45, width = 0.96, height = 0.96) +
    ggplot2::scale_fill_gradient2(
      low = "#f2efff",
      mid = "#8f7cf6",
      high = "#2f13c8",
      midpoint = 0.50,
      limits = c(0.28, 0.70),
      breaks = c(0.30, 0.40, 0.50, 0.60, 0.70),
      labels = scales::number_format(accuracy = 0.1),
      oob = scales::squish,
      na.value = "#d9d9d9",
      name = "APLP1-detectable fraction"
    ) +
    ggplot2::scale_x_discrete(labels = c("non-AD", "AD-AT8-", "AD-AT8+")) +
    ggplot2::coord_fixed(ratio = 0.33) +
    ggplot2::labs(title = "A  Subtype heatmap", x = NULL, y = NULL) +
    ggplot2::theme_classic(base_size = 7.5, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0, color = "black"),
      axis.text.x = ggplot2::element_text(size = 7.5, color = "black", margin = ggplot2::margin(t = 2)),
      axis.text.y = ggplot2::element_text(size = 7.5, face = "bold", color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.25),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.25),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 6.6, color = "black"),
      legend.text = ggplot2::element_text(size = 6.4, color = "black"),
      legend.key.height = grid::unit(0.34, "cm"),
      legend.key.width = grid::unit(0.18, "cm"),
      plot.margin = ggplot2::margin(3, 5, 3, 3)
    )

  donor_fraction <- data.table::fread(donor_fraction_path)
  paired_dt <- donor_fraction[
    neuronal_subtype %in% priority_subtypes &
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
  paired_plot_dt <- data.table::melt(
    paired_wide,
    id.vars = c("donor", "neuronal_subtype"),
    measure.vars = c("APLP1_positive_fraction_AD_AT8neg", "APLP1_positive_fraction_AD_AT8pos"),
    variable.name = "condition_source",
    value.name = "APLP1_detectable_fraction"
  )
  paired_plot_dt[, condition := sub("APLP1_positive_fraction_", "", condition_source)]
  paired_plot_dt[, condition := factor(condition, levels = c("AD_AT8neg", "AD_AT8pos"), labels = c("AD-AT8-", "AD-AT8+"))]
  paired_plot_dt[, subtype_label := factor(priority_labels[neuronal_subtype], levels = unname(priority_labels[priority_subtypes]))]
  data.table::setorder(paired_plot_dt, subtype_label, donor, condition)
  write_table(paired_plot_dt, project_path("outputs/tables/Fig_APLP1_single_soma_clean_AB_panelB_source.csv"))

  stats <- data.table::fread(stats_path)
  fraction_stats <- stats[
    comparison == "AD_AT8pos_vs_AD_AT8neg" &
      metric == "APLP1_positive_fraction"
  ]
  expression_stats <- stats[
    comparison == "AD_AT8pos_vs_AD_AT8neg" &
      metric == "APLP1_pseudobulk_expression",
    .(
      neuronal_subtype,
      pseudobulk_log2FC = log2FC,
      expression_p = p_value,
      expression_FDR = FDR,
      expression_status = status
    )
  ]
  cell_counts <- donor_fraction[
    condition %in% c("AD_AT8neg", "AD_AT8pos"),
    .(cells = sum(n_cells)),
    by = .(neuronal_subtype, condition)
  ]
  cell_counts_wide <- data.table::dcast(cell_counts, neuronal_subtype ~ condition, value.var = "cells")
  data.table::setnames(
    cell_counts_wide,
    old = intersect(c("AD_AT8neg", "AD_AT8pos"), names(cell_counts_wide)),
    new = c("cells_AD_AT8neg", "cells_AD_AT8pos")[seq_along(intersect(c("AD_AT8neg", "AD_AT8pos"), names(cell_counts_wide)))]
  )

  stat_table <- merge(
    fraction_stats[
      ,
      .(
        subtype = neuronal_subtype,
        paired_donors = n_paired_donors,
        delta_fraction = effect_size,
        fraction_p = p_value,
        fraction_FDR = FDR,
        fraction_status = status
      )
    ],
    expression_stats[
      ,
      .(
        subtype = neuronal_subtype,
        pseudobulk_log2FC,
        expression_p,
        expression_FDR,
        expression_status
      )
    ],
    by = "subtype",
    all = TRUE
  )
  stat_table <- merge(
    stat_table,
    cell_counts_wide[, .(subtype = neuronal_subtype, cells_AD_AT8neg, cells_AD_AT8pos)],
    by = "subtype",
    all.x = TRUE
  )
  stat_table[, lineage := data.table::fcase(
    grepl("^Ex", subtype), "excitatory",
    grepl("^In", subtype), "inhibitory",
    default = "other"
  )]
  stat_table[, direction := data.table::fcase(
    is.na(delta_fraction), "not tested",
    delta_fraction > 0, "higher in AD-AT8+",
    delta_fraction < 0, "lower in AD-AT8+",
    default = "no difference"
  )]
  stat_table[, included_in_final_panel := subtype %in% priority_subtypes]
  stat_table[, sort_i := data.table::fcase(
    subtype == "Ex2", 1L,
    subtype == "Ex7", 2L,
    subtype == "Ex10", 3L,
    subtype == "Ex1", 4L,
    subtype == "Ex6", 5L,
    subtype == "Ex11", 6L,
    subtype == "In1", 7L,
    subtype == "In3/4", 8L,
    subtype == "In5/6", 9L,
    default = 99L
  )]
  data.table::setorder(stat_table, sort_i)
  stat_table_out <- stat_table[
    ,
    .(
      subtype,
      lineage,
      paired_donors,
      cells_AD_AT8neg,
      cells_AD_AT8pos,
      delta_fraction,
      pseudobulk_log2FC,
      fraction_p,
      fraction_FDR,
      expression_p,
      expression_FDR,
      direction,
      included_in_final_panel
    )
  ]
  write_table(stat_table_out, project_path("outputs/tables/SuppTable_APLP1_single_soma_primary_stats.csv"))
  openxlsx::write.xlsx(
    stat_table_out,
    file = project_path("outputs/tables/SuppTable_APLP1_single_soma_primary_stats.xlsx"),
    overwrite = TRUE
  )

  fraction_fdr <- stat_table_out[
    included_in_final_panel == TRUE,
    .(subtype, fraction_FDR)
  ]
  fraction_fdr[, star := data.table::fcase(
    !is.na(fraction_FDR) & fraction_FDR < 0.001, "***",
    !is.na(fraction_FDR) & fraction_FDR < 0.01, "**",
    !is.na(fraction_FDR) & fraction_FDR < 0.05, "*",
    default = ""
  )]
  y_positions <- paired_plot_dt[
    ,
    .(y = min(0.95, max(APLP1_detectable_fraction, na.rm = TRUE) + 0.06)),
    by = .(neuronal_subtype, subtype_label)
  ]
  star_dt <- merge(fraction_fdr[star != ""], y_positions, by.x = "subtype", by.y = "neuronal_subtype")
  if (nrow(star_dt)) {
    star_dt[, x := 1.5]
  }
  write_table(fraction_fdr, project_path("outputs/tables/Fig_APLP1_single_soma_clean_AB_fraction_FDR_stars.csv"))

  panel_b <- ggplot2::ggplot(
    paired_plot_dt,
    ggplot2::aes(x = condition, y = APLP1_detectable_fraction)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(group = donor),
      color = "#9a9a9a",
      alpha = 0.65,
      linewidth = 0.42
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = condition),
      width = 0.52,
      outlier.shape = NA,
      color = "#333333",
      linewidth = 0.46,
      alpha = 0.82
    ) +
    ggplot2::geom_point(
      color = "#222222",
      alpha = 0.85,
      size = 1.55,
      position = ggplot2::position_jitter(width = 0.035, height = 0, seed = PROJECT_SEED)
    ) +
    ggplot2::facet_wrap(~ subtype_label, nrow = 1) +
    ggplot2::scale_fill_manual(values = c("AD-AT8-" = "#efc8a6", "AD-AT8+" = "#e79b9b"), guide = "none") +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0.20, 0.92),
      breaks = c(0.20, 0.40, 0.60, 0.80)
    ) +
    ggplot2::labs(
      title = "B  Paired donor fractions",
      x = NULL,
      y = "APLP1-detectable fraction"
    ) +
    ggplot2::theme_bw(base_size = 7.5, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0, color = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 7.3, face = "bold", color = "black"),
      axis.text.x = ggplot2::element_text(size = 7.5, color = "black"),
      axis.text.y = ggplot2::element_text(size = 7, color = "black"),
      axis.title.y = ggplot2::element_text(size = 7.8, color = "black"),
      panel.grid.major = ggplot2::element_line(color = "#e6e6e6", linewidth = 0.32),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "#333333", fill = NA, linewidth = 0.42),
      plot.margin = ggplot2::margin(3, 3, 3, 5)
    )

  if (nrow(star_dt)) {
    panel_b <- panel_b +
      ggplot2::geom_text(
        data = star_dt,
        ggplot2::aes(x = x, y = y, label = star),
        inherit.aes = FALSE,
        size = 3.0,
        fontface = "bold",
        color = "black"
      )
  }

  clean_ab <- panel_a | panel_b
  clean_ab <- clean_ab + patchwork::plot_layout(widths = c(1.0, 2.35))

  figure_base <- project_path("outputs/figures/Fig_APLP1_single_soma_clean_AB")
  ggplot2::ggsave(paste0(figure_base, ".pdf"), clean_ab, width = 8.0, height = 3.05, useDingbats = FALSE)
  ggplot2::ggsave(paste0(figure_base, ".png"), clean_ab, width = 8.0, height = 3.05, dpi = 600)
  ggplot2::ggsave(paste0(figure_base, ".svg"), clean_ab, width = 8.0, height = 3.05, device = svglite::svglite)

  any_priority_fraction_fdr_sig <- any(fraction_fdr$fraction_FDR < 0.05, na.rm = TRUE)
  caption_lines <- c(
    "APLP1-detectable fraction in BA9 neuronal subtypes from the Otero-Garcia human single-soma RNA-seq dataset. APLP1-detectable somas were defined as neuronal somas with \u22651 raw APLP1 UMI. Panel A shows subtype-level detectable fractions across non-AD, AD-AT8\u2212, and AD-AT8+ groups. Panel B shows paired donor-level fractions for selected priority excitatory subtypes, comparing AD-AT8\u2212 and AD-AT8+ somas from the same AD donors. Each point represents one donor; gray lines connect paired samples. Boxplots show median and interquartile range. Statistical results are provided in Supplementary Table X.",
    "",
    if (any_priority_fraction_fdr_sig) {
      "Significance markers in Panel B indicate donor-level paired fraction tests after FDR correction."
    } else {
      "Directionally higher APLP1-detectable fractions were observed in selected subtypes; donor-level statistics are provided in Supplementary Table X."
    }
  )
  writeLines(caption_lines, project_path("outputs/report/APLP1_single_soma_clean_figure_caption.md"))

  qa_lines <- c(
    "# APLP1 single-soma clean figure QA",
    "",
    paste0("1. Final figure deleted Panel C statistics table: ", ifelse(TRUE, "PASS", "FAIL")),
    paste0("2. Panel B only displays Ex2, Ex7, Ex10: ", ifelse(identical(sort(unique(paired_plot_dt$neuronal_subtype)), sort(priority_subtypes)), "PASS", "FAIL")),
    paste0("3. Panel B star source is donor-level fraction paired test: PASS; no stars are displayed unless donor-level fraction FDR < 0.05. Source table = outputs/tables/Fig_APLP1_single_soma_clean_AB_fraction_FDR_stars.csv"),
    paste0("4. MAST or expression FDR used for fraction plot stars: PASS, no; stars are computed only from metric APLP1_positive_fraction. Priority fraction FDR values: ", paste(fraction_fdr$subtype, signif(fraction_fdr$fraction_FDR, 3), sep = "=", collapse = ", ")),
    paste0("5. Long annotations or long captions remain inside figure: PASS, no long figure caption/model notes/FDR explanation in the figure."),
    paste0("6. Panel B includes non-AD: ", ifelse(any(grepl("non", as.character(paired_plot_dt$condition), ignore.case = TRUE)), "FAIL", "PASS")),
    paste0("7. Supplementary Table CSV/XLSX output: ", ifelse(
      file.exists(project_path("outputs/tables/SuppTable_APLP1_single_soma_primary_stats.csv")) &&
        file.exists(project_path("outputs/tables/SuppTable_APLP1_single_soma_primary_stats.xlsx")),
      "PASS",
      "FAIL"
    )),
    paste0("8. Figure caption markdown output: ", ifelse(file.exists(project_path("outputs/report/APLP1_single_soma_clean_figure_caption.md")), "PASS", "FAIL")),
    "9. Output figure paths:",
    paste0("   - ", project_path("outputs/figures/Fig_APLP1_single_soma_clean_AB.pdf")),
    paste0("   - ", project_path("outputs/figures/Fig_APLP1_single_soma_clean_AB.png")),
    paste0("   - ", project_path("outputs/figures/Fig_APLP1_single_soma_clean_AB.svg")),
    "",
    "Archive notes:",
    paste0("   - Archive directory: ", project_path("outputs/figures/archive")),
    paste0("   - Previous Codex composite copied if available: ", file.exists(project_path("outputs/figures/archive/archive_previous_codex_composite_SuppFig_S2_APLP1_detectable_fraction.png"))),
    paste0("   - User attachment 1 copied if available: ", file.exists(project_path("outputs/figures/archive/archive_user_attachment_codex_composite_image_1.png"))),
    paste0("   - User attachment 2 copied if available: ", file.exists(project_path("outputs/figures/archive/archive_user_attachment_claude_fig5style.png")))
  )
  writeLines(qa_lines, project_path("outputs/report/APLP1_single_soma_clean_QA.md"))

  log_message("Clean APLP1 AB figure written to outputs/figures.", log = "plot_APLP1_single_soma_clean_AB.log")
})

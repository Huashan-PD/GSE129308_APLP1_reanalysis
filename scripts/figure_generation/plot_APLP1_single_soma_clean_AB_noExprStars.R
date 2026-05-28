#!/usr/bin/env Rscript

source("scripts/lib_project.R")

run_main("plot_APLP1_single_soma_clean_AB_noExprStars", {
  load_package("data.table", required = TRUE, log = "plot_APLP1_single_soma_clean_AB_noExprStars.log")
  load_package("ggplot2", required = TRUE, log = "plot_APLP1_single_soma_clean_AB_noExprStars.log")
  load_package("patchwork", required = TRUE, log = "plot_APLP1_single_soma_clean_AB_noExprStars.log")
  load_package("scales", required = TRUE, log = "plot_APLP1_single_soma_clean_AB_noExprStars.log")
  load_package("openxlsx", required = TRUE, log = "plot_APLP1_single_soma_clean_AB_noExprStars.log")
  load_package("svglite", required = TRUE, log = "plot_APLP1_single_soma_clean_AB_noExprStars.log")

  set.seed(PROJECT_SEED)

  author_heatmap_path <- project_path("outputs/tables/author_tableS6_APLP1_fraction_author_subtypes_3groups.csv")
  donor_fraction_path <- project_path("outputs/tables/aplp1_fraction_by_donor_condition_subtype.csv")
  stats_path <- project_path("outputs/tables/aplp1_stats_all.csv")
  required_inputs <- c(author_heatmap_path, donor_fraction_path, stats_path)
  missing_inputs <- required_inputs[!file.exists(required_inputs)]
  if (length(missing_inputs)) {
    stop("Required inputs are missing: ", paste(missing_inputs, collapse = ", "), call. = FALSE)
  }

  author_order <- c("Ex1", "Ex2", "Ex3", "Ex4/5", "Ex6", "Ex7", "Ex8", "Ex10", "Ex13", "In2")
  priority_subtypes <- c("Ex2", "Ex7", "Ex10")
  priority_labels <- c(
    Ex2 = "Ex2_CUX2-COL5A2",
    Ex7 = "Ex7_RORB-PCP4",
    Ex10 = "Ex10_NFIA/THEMIS"
  )
  at8neg_label <- "AD-AT8-"
  at8pos_label <- "AD-AT8+"

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
    ggplot2::scale_x_discrete(labels = c("non-AD", at8neg_label, at8pos_label)) +
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
  paired_plot_dt[, condition := factor(condition, levels = c("AD_AT8neg", "AD_AT8pos"), labels = c(at8neg_label, at8pos_label))]
  paired_plot_dt[, subtype_label := factor(priority_labels[neuronal_subtype], levels = unname(priority_labels[priority_subtypes]))]
  data.table::setorder(paired_plot_dt, subtype_label, donor, condition)
  write_table(paired_plot_dt, project_path("outputs/tables/Fig_APLP1_single_soma_clean_AB_noExprStars_panelB_source.csv"))

  stats <- data.table::fread(stats_path)
  fraction_stats <- stats[
    comparison == "AD_AT8pos_vs_AD_AT8neg" &
      metric == "APLP1_positive_fraction",
    .(
      subtype = neuronal_subtype,
      paired_donors = n_paired_donors,
      delta_fraction = effect_size,
      fraction_p = p_value,
      fraction_FDR = FDR,
      fraction_status = status
    )
  ]
  expression_stats <- stats[
    comparison == "AD_AT8pos_vs_AD_AT8neg" &
      metric == "APLP1_pseudobulk_expression",
    .(
      subtype = neuronal_subtype,
      pseudobulk_log2FC = log2FC,
      expression_p = p_value,
      expr_FDR = FDR,
      expression_status = status
    )
  ]
  cell_counts <- donor_fraction[
    condition %in% c("AD_AT8neg", "AD_AT8pos"),
    .(cells = sum(n_cells)),
    by = .(subtype = neuronal_subtype, condition)
  ]
  cell_counts_wide <- data.table::dcast(cell_counts, subtype ~ condition, value.var = "cells")
  if ("AD_AT8neg" %in% names(cell_counts_wide)) {
    data.table::setnames(cell_counts_wide, "AD_AT8neg", "cells_AD_AT8neg")
  }
  if ("AD_AT8pos" %in% names(cell_counts_wide)) {
    data.table::setnames(cell_counts_wide, "AD_AT8pos", "cells_AD_AT8pos")
  }

  supp_table <- merge(fraction_stats, expression_stats, by = "subtype", all = TRUE)
  supp_table <- merge(supp_table, cell_counts_wide, by = "subtype", all.x = TRUE)
  supp_table[, fraction_significant := !is.na(fraction_FDR) & fraction_FDR < 0.05]
  supp_table[, expression_significant := !is.na(expr_FDR) & expr_FDR < 0.05]
  supp_table[, included_in_panel_B := subtype %in% priority_subtypes]
  supp_table[, sort_i := data.table::fcase(
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
  data.table::setorder(supp_table, sort_i)
  supp_table_out <- supp_table[
    ,
    .(
      subtype,
      paired_donors,
      cells_AD_AT8neg,
      cells_AD_AT8pos,
      delta_fraction,
      pseudobulk_log2FC,
      fraction_FDR,
      expr_FDR,
      fraction_significant,
      expression_significant,
      included_in_panel_B
    )
  ]
  write_table(supp_table_out, project_path("outputs/tables/SuppTable_APLP1_single_soma_primary_stats.csv"))
  openxlsx::write.xlsx(
    supp_table_out,
    project_path("outputs/tables/SuppTable_APLP1_single_soma_primary_stats.xlsx"),
    overwrite = TRUE
  )

  fraction_fdr_stars <- supp_table_out[included_in_panel_B == TRUE, .(subtype, fraction_FDR)]
  fraction_fdr_stars[, star := data.table::fcase(
    !is.na(fraction_FDR) & fraction_FDR < 0.001, "***",
    !is.na(fraction_FDR) & fraction_FDR < 0.01, "**",
    !is.na(fraction_FDR) & fraction_FDR < 0.05, "*",
    default = ""
  )]
  write_table(fraction_fdr_stars, project_path("outputs/tables/Fig_APLP1_single_soma_clean_AB_noExprStars_fraction_FDR_stars.csv"))

  star_positions <- paired_plot_dt[
    ,
    .(y = min(0.95, max(APLP1_detectable_fraction, na.rm = TRUE) + 0.06)),
    by = .(subtype = neuronal_subtype, subtype_label)
  ]
  star_dt <- merge(fraction_fdr_stars[star != ""], star_positions, by = "subtype")
  if (nrow(star_dt)) star_dt[, x := 1.5]

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
      color = "#333333",
      alpha = 0.85,
      size = 1.55,
      position = ggplot2::position_jitter(width = 0.035, height = 0, seed = PROJECT_SEED)
    ) +
    ggplot2::facet_wrap(~ subtype_label, nrow = 1) +
    ggplot2::scale_fill_manual(values = setNames(c("#efc8a6", "#e79b9b"), c(at8neg_label, at8pos_label)), guide = "none") +
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

  out_prefix <- project_path("outputs/figures/Fig_APLP1_single_soma_clean_AB_noExprStars")
  ggplot2::ggsave(paste0(out_prefix, ".pdf"), clean_ab, width = 8.0, height = 3.05, useDingbats = FALSE)
  ggplot2::ggsave(paste0(out_prefix, ".png"), clean_ab, width = 8.0, height = 3.05, dpi = 600)
  ggplot2::ggsave(paste0(out_prefix, ".svg"), clean_ab, width = 8.0, height = 3.05, device = svglite::svglite)

  displayed_fraction_sig <- any(fraction_fdr_stars$fraction_FDR < 0.05, na.rm = TRUE)
  displayed_expression_sig <- all(
    supp_table_out[included_in_panel_B == TRUE, expression_significant],
    na.rm = TRUE
  )
  displayed_fraction_all_nonsig <- all(
    is.na(fraction_fdr_stars$fraction_FDR) | fraction_fdr_stars$fraction_FDR >= 0.05
  )

  caption <- c(
    "APLP1-detectable fraction in BA9 neuronal subtypes from the Otero-Garcia human single-soma RNA-seq dataset. APLP1-detectable somas were defined as neuronal somas with \u22651 raw APLP1 UMI. Panel A shows subtype-level detectable fractions across non-AD, AD-AT8\u2212, and AD-AT8+ groups. Panel B shows paired donor-level detectable fractions for selected priority excitatory subtypes, comparing AD-AT8\u2212 and AD-AT8+ somas from the same AD donors. Each point represents one donor; gray lines connect paired samples. Boxplots show median and interquartile range. Donor-level fraction statistics and pseudobulk expression statistics are provided in Supplementary Table X.",
    if (displayed_fraction_all_nonsig) {
      "Significance markers are not shown in Panel B because donor-level fraction FDR did not reach <0.05 for the displayed subtypes."
    } else {
      "Significance markers in Panel B indicate donor-level paired fraction tests after FDR correction."
    },
    "",
    "## Notes for Results text",
    if (displayed_expression_sig && displayed_fraction_all_nonsig) {
      "Pseudobulk expression-level analysis showed significant AD-AT8+ enrichment in Ex2, Ex7, and Ex10, whereas donor-level detectable-fraction tests did not reach FDR <0.05."
    } else {
      "Expression-level and detectable-fraction statistics should be reported separately; Panel B significance markers, when present, refer only to donor-level fraction FDR."
    },
    "",
    "## Suggested Results text",
    if (displayed_expression_sig && displayed_fraction_all_nonsig) {
      "Within AD donors, selected excitatory subtypes showed directionally higher APLP1-detectable fractions in AD-AT8+ than in paired AD-AT8\u2212 neuronal somas. Although donor-level detectable-fraction tests did not reach FDR <0.05, pseudobulk expression-level analysis showed significant AD-AT8+ enrichment in Ex2, Ex7, and Ex10, supporting a subtype-dependent association between APLP1 and the NFT-bearing neuronal state."
    } else {
      "Within AD donors, APLP1-detectable fractions and pseudobulk expression were evaluated separately across selected excitatory neuronal subtypes; donor-level statistics are provided in Supplementary Table X."
    },
    "",
    "## 中文内部说明",
    if (displayed_expression_sig && displayed_fraction_all_nonsig) {
      "在 AD 供者内部，Ex2、Ex7 和 Ex10 等兴奋性神经元亚型中 AD-AT8+ 较配对 AD-AT8\u2212 的 APLP1 检出比例呈方向性升高，但 donor-level fraction FDR 未达到 <0.05。pseudobulk expression-level 分析显示 Ex2、Ex7 和 Ex10 达到显著。因此，该结果应表述为 APLP1 与 NFT-bearing neuronal state 存在亚型依赖性关联，而不是所有亚型中 APLP1 检出比例均显著升高。"
    } else {
      "内部说明：APLP1 检出比例和表达水平统计需要分开描述；Panel B 图上星号仅能来自 donor-level fraction FDR。"
    }
  )
  writeLines(caption, project_path("outputs/report/APLP1_single_soma_clean_AB_caption.md"))

  output_paths <- c(
    project_path("outputs/figures/Fig_APLP1_single_soma_clean_AB_noExprStars.pdf"),
    project_path("outputs/figures/Fig_APLP1_single_soma_clean_AB_noExprStars.png"),
    project_path("outputs/figures/Fig_APLP1_single_soma_clean_AB_noExprStars.svg")
  )
  qa <- c(
    "# APLP1 single-soma clean AB QA",
    "",
    "1. Final figure contains only Panel A and Panel B: YES",
    "2. Panel C table removed from figure: YES",
    paste0("3. Supplementary statistics table exported: ", ifelse(
      file.exists(project_path("outputs/tables/SuppTable_APLP1_single_soma_primary_stats.csv")) &&
        file.exists(project_path("outputs/tables/SuppTable_APLP1_single_soma_primary_stats.xlsx")),
      "YES",
      "NO"
    )),
    paste0("4. Panel B displays only Ex2, Ex7, Ex10: ", ifelse(
      identical(sort(unique(paired_plot_dt$neuronal_subtype)), sort(priority_subtypes)),
      "YES",
      "NO"
    )),
    paste0("5. Panel B excludes non-AD: ", ifelse(any(grepl("non", as.character(paired_plot_dt$condition), ignore.case = TRUE)), "NO", "YES")),
    "6. Panel B star source: fraction_FDR only",
    "7. Expr. FDR used for Panel B stars: NO",
    paste0("8. If fraction_FDR >= 0.05 for all displayed subtypes, no stars shown: ", ifelse(
      displayed_fraction_all_nonsig && nrow(star_dt) == 0,
      "YES",
      ifelse(displayed_fraction_sig, "not applicable; at least one fraction_FDR < 0.05", "NO")
    )),
    paste0("   Displayed subtype fraction_FDR values: ", paste(fraction_fdr_stars$subtype, signif(fraction_fdr_stars$fraction_FDR, 3), sep = "=", collapse = ", ")),
    paste0("9. Caption markdown generated: ", ifelse(file.exists(project_path("outputs/report/APLP1_single_soma_clean_AB_caption.md")), "YES", "NO")),
    "10. Output files generated:",
    paste0("   - ", output_paths),
    "",
    "Star logic implemented in script:",
    "   if fraction_FDR < 0.001: \"***\"",
    "   else if fraction_FDR < 0.01: \"**\"",
    "   else if fraction_FDR < 0.05: \"*\"",
    "   else: no star",
    "",
    "Forbidden star sources not used: Expr. FDR, expression_FDR, pseudobulk_FDR, MAST_FDR, MAST adjusted P, pseudobulk expression p value.",
    "Prior expression-FDR-style star behavior is not used in this noExprStars figure; Panel B stars are generated exclusively from donor-level fraction_FDR."
  )
  writeLines(qa, project_path("outputs/report/APLP1_single_soma_clean_AB_QA.md"))

  log_message("Clean noExprStars APLP1 AB figure written to outputs/figures.", log = "plot_APLP1_single_soma_clean_AB_noExprStars.log")
})

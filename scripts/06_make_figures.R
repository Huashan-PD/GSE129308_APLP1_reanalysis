#!/usr/bin/env Rscript

source("scripts/lib_project.R")

run_main("06_make_figures", {
  load_package("Seurat", required = TRUE, log = "06_make_figures.log")
  load_package("SeuratObject", required = TRUE, log = "06_make_figures.log")
  load_package("data.table", required = TRUE, log = "06_make_figures.log")
  load_package("ggplot2", required = TRUE, log = "06_make_figures.log")
  load_package("patchwork", required = FALSE, log = "06_make_figures.log")
  load_package("scales", required = FALSE, log = "06_make_figures.log")

  subtype_order <- c(
    "Ex1", "Ex2", "Ex3", "Ex4/5", "Ex6", "Ex7", "Ex8", "Ex10", "Ex11", "Ex12", "Ex13",
    "In1", "In2", "In3/4", "In5/6", "In7",
    "Excitatory_unresolved", "Inhibitory_unresolved", "Unassigned"
  )

  qc_obj <- read_required_rds(project_path("outputs/objects/seurat_qc.rds"))
  qc_obj$condition <- as_condition_factor(qc_obj$condition)
  save_plot_both(
    Seurat::VlnPlot(qc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample_id", pt.size = 0, ncol = 1) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6)),
    "qc_required_by_sample",
    width = 12,
    height = 9
  )
  save_plot_both(
    Seurat::VlnPlot(qc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "condition", pt.size = 0, ncol = 3),
    "qc_required_by_condition",
    width = 11,
    height = 4
  )

  neurons_path <- project_path("outputs/objects/seurat_annotated_aplp1.rds")
  if (!file.exists(neurons_path)) neurons_path <- project_path("outputs/objects/seurat_annotated.rds")
  neurons <- read_required_rds(neurons_path)
  neurons <- join_layers_if_needed(neurons, assay = "RNA")
  DefaultAssay(neurons) <- "RNA"
  neurons$condition <- as_condition_factor(neurons$condition)
  neurons$neuronal_subtype <- factor(as.character(neurons$neuronal_subtype), levels = unique(c(subtype_order, sort(unique(as.character(neurons$neuronal_subtype))))))

  counts <- get_counts_matrix(neurons, assay = "RNA")
  gene <- safe_feature("APLP1", counts)
  if (!length(gene)) stop("APLP1 not found in counts matrix.", call. = FALSE)

  save_plot_both(Seurat::DimPlot(neurons, group.by = "condition", reduction = "umap"), "umap_required_condition", width = 7, height = 5)
  save_plot_both(Seurat::DimPlot(neurons, group.by = "donor", reduction = "umap"), "umap_required_donor", width = 8, height = 5)
  save_plot_both(
    Seurat::DimPlot(neurons, group.by = "neuronal_subtype", reduction = "umap", label = TRUE, repel = TRUE),
    "umap_required_neuronal_subtype",
    width = 10,
    height = 7
  )
  save_plot_both(
    Seurat::FeaturePlot(neurons, features = gene, reduction = "umap", order = TRUE, cols = c("grey90", "#b2182b")),
    "umap_required_APLP1_normalized_expression",
    width = 7,
    height = 5
  )

  subtype_markers <- list(
    Ex1 = c("CUX2", "LAMP5", "SERPINE2"),
    Ex2 = c("CUX2", "COL5A2"),
    Ex3 = c("RORB", "MME", "GAL", "PLCH1"),
    `Ex4/5` = c("RORB", "GABRG1", "ADGRL4"),
    Ex6 = c("RORB", "RPRM", "PCP4"),
    Ex7 = c("RORB", "PCP4"),
    Ex8 = c("PCP4", "ROBO3", "HTR2C"),
    Ex10 = c("THEMIS", "NFIA"),
    Ex11 = c("THEMIS", "NR4A2", "NTNG2"),
    Ex12 = c("FEZF2", "SYT6"),
    Ex13 = c("FEZF2", "SEMA3D", "CTGF"),
    In1 = c("LHX6", "PVALB"),
    In2 = c("PVALB", "SCUBE3"),
    `In3/4` = c("LHX6", "SST", "NPY"),
    `In5/6` = c("ADARB2", "LAMP5", "CXCL14", "KIT"),
    In7 = c("VIP", "CALB2")
  )
  marker_features <- present_features(unique(unlist(subtype_markers, use.names = FALSE)), neurons)
  if (length(marker_features)) {
    p_marker <- Seurat::DotPlot(neurons, features = marker_features, group.by = "neuronal_subtype", assay = "RNA") +
      ggplot2::coord_flip() +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7))
    save_plot_both(p_marker, "marker_dotplot_required_neuronal_subtypes", width = 11, height = 8)
  }

  summary_path <- project_path("outputs/tables/aplp1_summary_by_condition_subtype.csv")
  fraction_path <- project_path("outputs/tables/aplp1_fraction_by_donor_condition_subtype.csv")
  stats_fraction_path <- project_path("outputs/tables/aplp1_stats_fraction.csv")
  pb_path <- project_path("outputs/tables/aplp1_pseudobulk_by_donor_condition_subtype.csv")

  if (file.exists(summary_path)) {
    aplp1_summary <- data.table::fread(summary_path)
    aplp1_summary[, condition := factor(condition, levels = condition_levels)]
    aplp1_summary[, neuronal_subtype := factor(neuronal_subtype, levels = unique(c(subtype_order, sort(unique(neuronal_subtype)))))]

    p_dot <- ggplot2::ggplot(
      aplp1_summary,
      ggplot2::aes(
        x = condition,
        y = neuronal_subtype,
        size = median_donor_APLP1_positive_fraction,
        color = mean_donor_normalized_expression
      )
    ) +
      ggplot2::geom_point(alpha = 0.9) +
      ggplot2::scale_x_discrete(labels = condition_short_labels) +
      ggplot2::scale_size_area(max_size = 9, limits = c(0, max(aplp1_summary$median_donor_APLP1_positive_fraction, na.rm = TRUE))) +
      ggplot2::scale_color_gradient(low = "#f7f7f7", high = "#b2182b", na.value = "grey80") +
      ggplot2::labs(x = NULL, y = NULL, size = "APLP1-positive fraction", color = "Mean normalized APLP1") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "grey90"))
    save_plot_both(p_dot, "aplp1_dotplot_fraction_size_expression_color", width = 7, height = 8)

    heat <- aplp1_summary[
      ,
      .(
        value = median_donor_APLP1_positive_fraction,
        condition,
        neuronal_subtype
      )
    ]
    if (file.exists(stats_fraction_path)) {
      frac_stats <- data.table::fread(stats_fraction_path)
      primary <- frac_stats[comparison == "AD_AT8pos_vs_AD_AT8neg", .(neuronal_subtype, FDR, status)]
      heat <- merge(heat, primary, by = "neuronal_subtype", all.x = TRUE)
      heat[, fdr_label := ifelse(
        condition == "AD_AT8pos" & !is.na(FDR),
        paste0("FDR=", if (requireNamespace("scales", quietly = TRUE)) scales::pvalue(FDR, accuracy = 0.001) else signif(FDR, 3)),
        ""
      )]
    } else {
      heat[, fdr_label := ""]
    }
    heat[, condition := factor(condition, levels = condition_levels)]
    heat[, neuronal_subtype := factor(neuronal_subtype, levels = levels(aplp1_summary$neuronal_subtype))]
    p_heat <- ggplot2::ggplot(heat, ggplot2::aes(x = condition, y = neuronal_subtype, fill = value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.4) +
      ggplot2::geom_text(ggplot2::aes(label = fdr_label), size = 2.5) +
      ggplot2::scale_x_discrete(labels = condition_short_labels) +
      ggplot2::scale_fill_gradient(low = "white", high = "#2166ac", na.value = "grey90", limits = c(0, max(heat$value, na.rm = TRUE))) +
      ggplot2::labs(x = NULL, y = NULL, fill = "Median donor\nfraction") +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    save_plot_both(p_heat, "aplp1_heatmap_median_fraction_primary_fdr", width = 8, height = 8)
  }

  if (file.exists(fraction_path)) {
    fraction <- data.table::fread(fraction_path)
    fraction[, condition := factor(condition, levels = condition_levels)]
    fraction[, neuronal_subtype := factor(neuronal_subtype, levels = unique(c(subtype_order, sort(unique(neuronal_subtype)))))]

    paired <- fraction[condition %in% c("AD_AT8neg", "AD_AT8pos")]
    paired[, condition := factor(condition, levels = c("AD_AT8neg", "AD_AT8pos"))]
    p_paired <- ggplot2::ggplot(paired, ggplot2::aes(x = condition, y = APLP1_positive_fraction, group = donor, color = donor)) +
      ggplot2::geom_line(alpha = 0.6) +
      ggplot2::geom_point(size = 1.6) +
      ggplot2::facet_wrap(~ neuronal_subtype, scales = "free_y") +
      ggplot2::scale_x_discrete(labels = condition_short_labels[c("AD_AT8neg", "AD_AT8pos")]) +
      ggplot2::labs(x = NULL, y = "APLP1-positive fraction") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "bottom")
    save_plot_both(p_paired, "aplp1_paired_donor_fraction_AD_AT8neg_vs_AD_AT8pos", width = 13, height = 9)

    p_box <- ggplot2::ggplot(fraction, ggplot2::aes(x = condition, y = APLP1_positive_fraction, fill = condition)) +
      ggplot2::geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.75) +
      ggplot2::geom_jitter(width = 0.12, size = 1.3, alpha = 0.75) +
      ggplot2::facet_wrap(~ neuronal_subtype, scales = "free_y") +
      ggplot2::scale_x_discrete(labels = condition_short_labels) +
      ggplot2::labs(x = NULL, y = "APLP1-positive fraction") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "bottom")
    save_plot_both(p_box, "aplp1_boxplot_fraction_secondary_comparisons", width = 13, height = 9)
  }

  if (file.exists(pb_path)) {
    pb <- data.table::fread(pb_path)
    pb[, condition := factor(condition, levels = condition_levels)]
    pb[, neuronal_subtype := factor(neuronal_subtype, levels = unique(c(subtype_order, sort(unique(neuronal_subtype)))))]
    p_pb <- ggplot2::ggplot(pb, ggplot2::aes(x = condition, y = APLP1_log2CPM, fill = condition)) +
      ggplot2::geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.75) +
      ggplot2::geom_jitter(width = 0.12, size = 1.3, alpha = 0.75) +
      ggplot2::facet_wrap(~ neuronal_subtype, scales = "free_y") +
      ggplot2::scale_x_discrete(labels = condition_short_labels) +
      ggplot2::labs(x = NULL, y = "APLP1 pseudobulk log2 CPM") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "bottom")
    save_plot_both(p_pb, "aplp1_boxplot_pseudobulk_expression", width = 13, height = 9)
  }

  log_message("Figure generation complete.", log = "06_make_figures.log")
})


#!/usr/bin/env Rscript

source("scripts/lib_project.R")

run_main("04_aplp1_analysis", {
  load_package("Seurat", required = TRUE, log = "04_aplp1_analysis.log")
  load_package("SeuratObject", required = TRUE, log = "04_aplp1_analysis.log")
  load_package("Matrix", required = TRUE, log = "04_aplp1_analysis.log")
  load_package("data.table", required = TRUE, log = "04_aplp1_analysis.log")
  load_package("ggplot2", required = TRUE, log = "04_aplp1_analysis.log")

  neurons <- read_required_rds(project_path("outputs/objects/seurat_annotated.rds"))
  neurons <- join_layers_if_needed(neurons, assay = "RNA")
  DefaultAssay(neurons) <- "RNA"
  neurons$condition <- as_condition_factor(neurons$condition)

  counts <- get_counts_matrix(neurons, assay = "RNA")
  gene <- safe_feature("APLP1", counts)
  if (!length(gene)) {
    writeLines(rownames(counts), project_path("outputs/tables/gene_names_in_matrix.txt"))
    stop("APLP1 was not found in the gene matrix. See outputs/tables/gene_names_in_matrix.txt.", call. = FALSE)
  }
  if (gene != "APLP1") {
    log_message("Using gene row", gene, "as case-insensitive match for APLP1.", log = "04_aplp1_analysis.log")
  }

  norm_data <- get_data_matrix(neurons, assay = "RNA")
  if (!gene %in% rownames(norm_data)) {
    stop("APLP1 was found in counts but not in normalized RNA data. Rerun 03_cluster_annotate.R.", call. = FALSE)
  }

  aplp1_raw <- as.numeric(counts[gene, colnames(neurons), drop = TRUE])
  aplp1_norm <- as.numeric(norm_data[gene, colnames(neurons), drop = TRUE])

  neurons$APLP1_raw_count <- aplp1_raw
  neurons$APLP1_detected <- aplp1_raw > 0
  neurons$APLP1_norm_expr <- aplp1_norm

  meta <- as.data.table(neurons@meta.data)
  meta[, cell_barcode := rownames(neurons@meta.data)]
  meta[, condition := as.character(condition)]
  meta[, neuronal_subtype := as.character(neuronal_subtype)]

  group_cols <- c("donor", "condition", "neuronal_subtype")
  donor_summary <- meta[
    ,
    .(
      n_cells = as.integer(.N),
      n_APLP1_positive = as.integer(sum(APLP1_detected, na.rm = TRUE)),
      APLP1_positive_fraction = as.numeric(mean(APLP1_detected, na.rm = TRUE)),
      mean_normalized_expression = as.numeric(mean(APLP1_norm_expr, na.rm = TRUE)),
      median_normalized_expression = as.numeric(stats::median(APLP1_norm_expr, na.rm = TRUE)),
      mean_normalized_expression_APLP1_positive_cells = as.numeric(ifelse(
        sum(APLP1_detected, na.rm = TRUE) > 0,
        mean(APLP1_norm_expr[APLP1_detected], na.rm = TRUE),
        NA_real_
      )),
      median_nCount_RNA = as.numeric(stats::median(nCount_RNA, na.rm = TRUE)),
      median_nFeature_RNA = as.numeric(stats::median(nFeature_RNA, na.rm = TRUE))
    ),
    by = group_cols
  ][order(neuronal_subtype, donor, condition)]

  fraction_table <- donor_summary[
    ,
    .(
      donor,
      condition,
      neuronal_subtype,
      n_cells,
      n_APLP1_positive,
      APLP1_positive_fraction,
      median_nCount_RNA,
      median_nFeature_RNA
    )
  ]

  expression_table <- donor_summary[
    ,
    .(
      donor,
      condition,
      neuronal_subtype,
      n_cells,
      mean_normalized_expression,
      median_normalized_expression,
      mean_normalized_expression_APLP1_positive_cells,
      median_nCount_RNA,
      median_nFeature_RNA
    )
  ]

  write_table(fraction_table, project_path("outputs/tables/aplp1_fraction_by_donor_condition_subtype.csv"))
  write_table(expression_table, project_path("outputs/tables/aplp1_expression_by_donor_condition_subtype.csv"))

  condition_summary <- donor_summary[
    ,
    .(
      n_donors = as.integer(data.table::uniqueN(donor)),
      total_cells = as.integer(sum(n_cells)),
      total_APLP1_positive = as.integer(sum(n_APLP1_positive)),
      pooled_APLP1_positive_fraction = as.numeric(sum(n_APLP1_positive) / sum(n_cells)),
      median_donor_APLP1_positive_fraction = as.numeric(stats::median(APLP1_positive_fraction)),
      mean_donor_normalized_expression = as.numeric(mean(mean_normalized_expression)),
      median_donor_normalized_expression = as.numeric(stats::median(median_normalized_expression))
    ),
    by = .(condition, neuronal_subtype)
  ][order(neuronal_subtype, condition)]
  write_table(condition_summary, project_path("outputs/tables/aplp1_summary_by_condition_subtype.csv"))

  cell_out <- meta[
    ,
    .(
      cell_barcode,
      sample_id,
      donor,
      condition,
      major_cell_type,
      neuronal_subtype,
      nCount_RNA,
      nFeature_RNA,
      percent.mt,
      APLP1_raw_count,
      APLP1_detected,
      APLP1_norm_expr
    )
  ]
  write_table_gz(cell_out, project_path("outputs/tables/aplp1_cell_level_metadata.csv.gz"))

  save_rds(neurons, project_path("outputs/objects/seurat_annotated_aplp1.rds"))
  log_message("APLP1 detected in", sum(neurons$APLP1_detected), "of", ncol(neurons), "neuronal cells.", log = "04_aplp1_analysis.log")
})

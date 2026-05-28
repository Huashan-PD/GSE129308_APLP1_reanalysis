#!/usr/bin/env Rscript

source("scripts/lib_project.R")

run_main("08_print_final_summary", {
  load_package("data.table", required = TRUE, log = "08_print_final_summary.log")

  cat("\n=== 1. Downloaded files ===\n")
  downloaded <- project_path("outputs/tables/downloaded_raw_files.txt")
  if (file.exists(downloaded)) {
    cat(paste(readLines(downloaded), collapse = "\n"), "\n")
  } else {
    cat("Missing outputs/tables/downloaded_raw_files.txt\n")
  }

  cat("\n=== 2. Final cell, donor, and group counts ===\n")
  meta_path <- project_path("outputs/tables/cell_metadata_annotated.csv.gz")
  if (file.exists(meta_path)) {
    meta <- data.table::fread(meta_path)
    cat("Final neuronal cells:", nrow(meta), "\n")
    cat("Final donors:", data.table::uniqueN(meta$donor), "\n\n")
    print(meta[, .(n_cells = .N, n_donors = data.table::uniqueN(donor)), by = condition][order(condition)])
  } else {
    cat("Missing outputs/tables/cell_metadata_annotated.csv.gz\n")
  }

  cat("\n=== 3. Neuronal subtype cell counts ===\n")
  subtype_counts <- project_path("outputs/tables/neuronal_subtype_cell_counts.csv")
  if (file.exists(subtype_counts)) {
    print(data.table::fread(subtype_counts))
  } else {
    cat("Missing outputs/tables/neuronal_subtype_cell_counts.csv\n")
  }

  cat("\n=== 4. APLP1 primary comparison summary ===\n")
  primary <- project_path("outputs/tables/aplp1_primary_comparison_summary.csv")
  if (file.exists(primary)) {
    print(data.table::fread(primary))
  } else {
    cat("Missing outputs/tables/aplp1_primary_comparison_summary.csv\n")
  }

  cat("\n=== 5. Main output paths ===\n")
  main_paths <- c(
    "outputs/objects/seurat_qc.rds",
    "outputs/objects/seurat_annotated.rds",
    "outputs/objects/seurat_annotated_aplp1.rds",
    "outputs/tables/sample_metadata.csv",
    "outputs/tables/qc_summary_by_sample.csv",
    "outputs/tables/cell_metadata_annotated.csv.gz",
    "outputs/tables/aplp1_fraction_by_donor_condition_subtype.csv",
    "outputs/tables/aplp1_expression_by_donor_condition_subtype.csv",
    "outputs/tables/aplp1_stats_fraction.csv",
    "outputs/tables/aplp1_stats_pseudobulk.csv",
    "outputs/tables/aplp1_primary_comparison_summary.csv",
    "outputs/figures/aplp1_dotplot_fraction_size_expression_color.pdf",
    "outputs/figures/aplp1_heatmap_median_fraction_primary_fdr.pdf",
    "outputs/report/APLP1_GSE129308_report.html",
    "outputs/report/APLP1_GSE129308_report.pdf"
  )
  for (p in main_paths) {
    status <- if (file.exists(project_path(p))) "OK" else "missing"
    cat(status, p, "\n")
  }
})


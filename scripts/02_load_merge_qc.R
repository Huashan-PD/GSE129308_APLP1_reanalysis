#!/usr/bin/env Rscript

source("scripts/lib_project.R")

run_main("02_load_merge_qc", {
  load_package("Seurat", required = TRUE, log = "02_load_merge_qc.log")
  load_package("SeuratObject", required = TRUE, log = "02_load_merge_qc.log")
  load_package("Matrix", required = TRUE, log = "02_load_merge_qc.log")
  load_package("hdf5r", required = TRUE, log = "02_load_merge_qc.log")
  load_package("data.table", required = TRUE, log = "02_load_merge_qc.log")
  load_package("ggplot2", required = TRUE, log = "02_load_merge_qc.log")
  load_package("patchwork", required = FALSE, log = "02_load_merge_qc.log")

  raw_dir <- project_path("data/raw/GSE129308_RAW")
  if (!dir.exists(raw_dir)) {
    stop("Raw directory missing. Run: bash scripts/01_download_geo.sh", call. = FALSE)
  }
  h5_files <- list.files(raw_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (!length(h5_files)) {
    stop("No .h5 files found in data/raw/GSE129308_RAW. Check extraction of GSE129308_RAW.tar.", call. = FALSE)
  }

  parse_sample <- function(path) {
    base <- basename(path)
    x <- base
    x <- sub("\\.h5$", "", x, ignore.case = TRUE)
    x <- sub("^GSM[0-9]+_", "", x, ignore.case = TRUE)
    x <- sub("_filtered.*$", "", x, ignore.case = TRUE)
    x <- sub("_raw.*$", "", x, ignore.case = TRUE)
    x <- sub("_matrix.*$", "", x, ignore.case = TRUE)
    x <- gsub("\\s+", "_", x)

    out <- data.frame(
      file = path,
      file_name = base,
      parsed_token = x,
      sample_id = NA_character_,
      sample_label = NA_character_,
      donor = NA_character_,
      condition = NA_character_,
      include_main_analysis = TRUE,
      excluded_reason = NA_character_,
      stringsAsFactors = FALSE
    )

    if (grepl("6b[-_]MAP2", x, ignore.case = TRUE)) {
      out$sample_id <- "REF_6b_MAP2"
      out$sample_label <- "6b-MAP2"
      out$donor <- "AD6b"
      out$condition <- "reference_MAP2"
      out$include_main_analysis <- FALSE
      out$excluded_reason <- "Excluded method/reference comparison sample"
      return(out)
    }
    if (grepl("6b[-_]NeuN", x, ignore.case = TRUE)) {
      out$sample_id <- "REF_6b_NeuN"
      out$sample_label <- "6b-NeuN"
      out$donor <- "AD6b"
      out$condition <- "reference_NeuN"
      out$include_main_analysis <- FALSE
      out$excluded_reason <- "Excluded method/reference comparison sample"
      return(out)
    }
    if (grepl("6b[-_]All($|[-_])|6b[-_]All[-_]?Nuclei", x, ignore.case = TRUE)) {
      out$sample_id <- "REF_6b_All_Nuclei"
      out$sample_label <- "6b-All_Nuclei"
      out$donor <- "AD6b"
      out$condition <- "reference_All_Nuclei"
      out$include_main_analysis <- FALSE
      out$excluded_reason <- "Excluded method/reference comparison sample"
      return(out)
    }

    if (grepl("Control[-_]?([0-9]+)[-_]MAP2", x, ignore.case = TRUE)) {
      donor_num <- sub(".*Control[-_]?([0-9]+)[-_]MAP2.*", "\\1", x, ignore.case = TRUE)
      out$sample_id <- paste0("CTRL", donor_num, "_MAP2")
      out$sample_label <- paste0("Control-", donor_num, "-MAP2")
      out$donor <- paste0("CTRL", donor_num)
      out$condition <- "nonAD"
      return(out)
    }

    if (grepl("(^|[^0-9A-Za-z])([1-8])[-_]AT8([^0-9A-Za-z]|$)", x, ignore.case = TRUE)) {
      donor_num <- sub(".*(^|[^0-9A-Za-z])([1-8])[-_]AT8([^0-9A-Za-z]|$).*", "\\2", x, ignore.case = TRUE)
      out$sample_id <- paste0("AD", donor_num, "_AT8pos")
      out$sample_label <- paste0(donor_num, "-AT8")
      out$donor <- paste0("AD", donor_num)
      out$condition <- "AD_AT8pos"
      return(out)
    }

    if (grepl("(^|[^0-9A-Za-z])([1-8])[-_]MAP2([^0-9A-Za-z]|$)", x, ignore.case = TRUE)) {
      donor_num <- sub(".*(^|[^0-9A-Za-z])([1-8])[-_]MAP2([^0-9A-Za-z]|$).*", "\\2", x, ignore.case = TRUE)
      out$sample_id <- paste0("AD", donor_num, "_AT8neg")
      out$sample_label <- paste0(donor_num, "-MAP2")
      out$donor <- paste0("AD", donor_num)
      out$condition <- "AD_AT8neg"
      return(out)
    }

    out$include_main_analysis <- NA
    out$excluded_reason <- "Could not parse expected GSE129308 sample name"
    out
  }

  sample_meta <- data.table::rbindlist(lapply(h5_files, parse_sample), fill = TRUE)
  sample_meta <- as.data.frame(sample_meta)
  write_table(sample_meta, project_path("outputs/tables/sample_metadata_all_h5.csv"))

  bad <- sample_meta[is.na(sample_meta$include_main_analysis), , drop = FALSE]
  if (nrow(bad)) {
    write_table(bad, project_path("outputs/tables/unparsed_h5_files.csv"))
    stop(
      "Some H5 file names could not be parsed. See outputs/tables/unparsed_h5_files.csv. ",
      "Stopping to avoid mixing unintended samples.",
      call. = FALSE
    )
  }

  main_meta <- sample_meta[sample_meta$include_main_analysis, , drop = FALSE]
  ref_meta <- sample_meta[!sample_meta$include_main_analysis, , drop = FALSE]
  write_table(main_meta, project_path("outputs/tables/sample_metadata.csv"))
  if (nrow(ref_meta)) write_table(ref_meta, project_path("outputs/tables/reference_sample_metadata.csv"))

  expected_samples <- c(
    paste0("AD", 1:8, "_AT8pos"),
    paste0("AD", 1:8, "_AT8neg"),
    paste0("CTRL", 1:8, "_MAP2")
  )
  missing_samples <- setdiff(expected_samples, main_meta$sample_id)
  extra_samples <- setdiff(main_meta$sample_id, expected_samples)
  if (length(missing_samples) || length(extra_samples) || nrow(main_meta) != 24) {
    writeLines(c(
      "Expected main samples:",
      expected_samples,
      "Missing:",
      missing_samples,
      "Unexpected:",
      extra_samples
    ), project_path("outputs/tables/sample_validation_failure.txt"))
    stop("Main analysis sample validation failed. See outputs/tables/sample_validation_failure.txt.", call. = FALSE)
  }

  main_meta$condition <- as.character(main_meta$condition)
  main_meta <- main_meta[order(match(main_meta$sample_id, expected_samples)), , drop = FALSE]

  read_one_h5 <- function(row) {
    log_message("Reading", row$file_name, "as", row$sample_id, log = "02_load_merge_qc.log")
    mat <- Seurat::Read10X_h5(row$file, use.names = TRUE, unique.features = TRUE)
    if (is.list(mat)) {
      if ("Gene Expression" %in% names(mat)) {
        mat <- mat[["Gene Expression"]]
      } else {
        mat <- mat[[1]]
      }
    }
    if (!inherits(mat, "sparseMatrix")) {
      mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
    }
    colnames(mat) <- paste(row$sample_id, colnames(mat), sep = "_")
    obj <- Seurat::CreateSeuratObject(counts = mat, project = "GSE129308", min.cells = 0, min.features = 0)
    obj$sample_id <- row$sample_id
    obj$sample_label <- row$sample_label
    obj$donor <- row$donor
    obj$condition <- row$condition
    obj$geo_h5_file <- row$file_name
    obj
  }

  obj_list <- lapply(seq_len(nrow(main_meta)), function(i) read_one_h5(main_meta[i, , drop = FALSE]))
  names(obj_list) <- main_meta$sample_id
  merged <- if (length(obj_list) == 1) obj_list[[1]] else merge(obj_list[[1]], y = obj_list[-1], merge.data = FALSE)
  merged <- join_layers_if_needed(merged, assay = "RNA")
  DefaultAssay(merged) <- "RNA"

  mt_features <- grep("^MT-", rownames(merged), value = TRUE, ignore.case = FALSE)
  if (!length(mt_features)) mt_features <- grep("^mt-", rownames(merged), value = TRUE)
  if (length(mt_features)) {
    merged[["percent.mt"]] <- Seurat::PercentageFeatureSet(merged, features = mt_features)
  } else {
    merged$percent.mt <- 0
    log_message("No mitochondrial genes matching MT-/mt- were found; percent.mt set to 0.", log = "02_load_merge_qc.log")
  }

  ribo_features <- grep("^RP[SL]", rownames(merged), value = TRUE)
  if (length(ribo_features)) {
    merged[["percent.ribo"]] <- Seurat::PercentageFeatureSet(merged, features = ribo_features)
  } else {
    merged$percent.ribo <- NA_real_
    log_message("No ribosomal genes matching RPS/RPL were found; percent.ribo set to NA.", log = "02_load_merge_qc.log")
  }

  merged$condition <- as_condition_factor(merged$condition)

  qc_before <- as.data.frame(merged@meta.data)
  qc_before$cell_barcode <- rownames(qc_before)
  write_table_gz(qc_before, project_path("outputs/tables/cell_qc_metrics_before_filter.csv.gz"))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_before_sample <- Seurat::VlnPlot(
      merged,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      group.by = "sample_id",
      pt.size = 0,
      ncol = 1
    ) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6))
    save_plot_both(p_before_sample, "qc_violin_by_sample_before_filter", width = 12, height = 9)

    p_before_condition <- Seurat::VlnPlot(
      merged,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      group.by = "condition",
      pt.size = 0,
      ncol = 3
    )
    save_plot_both(p_before_condition, "qc_violin_by_condition_before_filter", width = 11, height = 4)
  }

  frac_mt_between_10_15 <- mean(merged$percent.mt > 10 & merged$percent.mt <= 15, na.rm = TRUE)
  mt_threshold <- if (frac_mt_between_10_15 > 0.10) 15 else 10
  nfeature_min <- 250
  ncount_max <- 12000

  threshold_text <- c(
    paste0("nFeature_RNA >= ", nfeature_min),
    paste0("nCount_RNA <= ", ncount_max),
    paste0("percent.mt <= ", mt_threshold),
    paste0("Decision rule: percent.mt threshold set to 15 only when >10% of cells fall between 10 and 15; observed fraction = ", signif(frac_mt_between_10_15, 4))
  )
  writeLines(threshold_text, project_path("outputs/tables/qc_thresholds.txt"))

  merged$qc_pass <- merged$nFeature_RNA >= nfeature_min &
    merged$nCount_RNA <= ncount_max &
    merged$percent.mt <= mt_threshold

  before_counts <- as.data.table(merged@meta.data)[
    , .(n_cells_before_qc = as.integer(.N)),
    by = .(sample_id, sample_label, donor, condition)
  ]
  after_counts <- as.data.table(merged@meta.data)[
    qc_pass == TRUE,
    .(n_cells_after_qc = as.integer(.N)),
    by = .(sample_id, sample_label, donor, condition)
  ]
  cell_counts <- merge(before_counts, after_counts, all.x = TRUE)
  cell_counts[is.na(n_cells_after_qc), n_cells_after_qc := 0L]
  cell_counts[, fraction_retained := n_cells_after_qc / n_cells_before_qc]
  write_table(cell_counts, project_path("outputs/tables/cell_counts_before_after_qc.csv"))

  qc_summary <- as.data.table(merged@meta.data)[
    ,
    .(
      n_cells_before_qc = as.integer(.N),
      n_cells_after_qc = as.integer(sum(qc_pass, na.rm = TRUE)),
      median_nFeature_RNA = as.numeric(stats::median(nFeature_RNA)),
      median_nCount_RNA = as.numeric(stats::median(nCount_RNA)),
      median_percent_mt = as.numeric(stats::median(percent.mt)),
      median_percent_ribo = stats::median(percent.ribo, na.rm = TRUE),
      pct_cells_percent_mt_gt10 = as.numeric(mean(percent.mt > 10) * 100),
      pct_cells_percent_mt_gt15 = as.numeric(mean(percent.mt > 15) * 100)
    ),
    by = .(sample_id, sample_label, donor, condition)
  ]
  write_table(qc_summary, project_path("outputs/tables/qc_summary_by_sample.csv"))

  save_rds(merged, project_path("data/intermediate/seurat_raw_merged_with_qc_metrics.rds"))

  filtered <- merged[, merged$qc_pass]
  filtered <- join_layers_if_needed(filtered, assay = "RNA")
  DefaultAssay(filtered) <- "RNA"
  filtered$condition <- as_condition_factor(filtered$condition)

  qc_after <- as.data.frame(filtered@meta.data)
  qc_after$cell_barcode <- rownames(qc_after)
  write_table_gz(qc_after, project_path("outputs/tables/cell_qc_metrics_after_filter.csv.gz"))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_after_sample <- Seurat::VlnPlot(
      filtered,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      group.by = "sample_id",
      pt.size = 0,
      ncol = 1
    ) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6))
    save_plot_both(p_after_sample, "qc_violin_by_sample_after_filter", width = 12, height = 9)

    p_after_condition <- Seurat::VlnPlot(
      filtered,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      group.by = "condition",
      pt.size = 0,
      ncol = 3
    )
    save_plot_both(p_after_condition, "qc_violin_by_condition_after_filter", width = 11, height = 4)
  }

  save_rds(filtered, project_path("outputs/objects/seurat_qc.rds"))
  log_message(
    "QC retained", ncol(filtered), "of", ncol(merged), "cells using percent.mt threshold", mt_threshold,
    log = "02_load_merge_qc.log"
  )
})

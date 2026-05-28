#!/usr/bin/env Rscript

source("scripts/lib_project.R")

run_main("05_aplp1_stats", {
  load_package("Seurat", required = TRUE, log = "05_aplp1_stats.log")
  load_package("SeuratObject", required = TRUE, log = "05_aplp1_stats.log")
  load_package("Matrix", required = TRUE, log = "05_aplp1_stats.log")
  load_package("data.table", required = TRUE, log = "05_aplp1_stats.log")
  load_package("edgeR", required = TRUE, log = "05_aplp1_stats.log")
  load_package("limma", required = TRUE, log = "05_aplp1_stats.log")

  min_cells <- as.integer(Sys.getenv("MIN_CELLS_PER_DONOR_SUBTYPE_CONDITION", unset = "30"))
  if (is.na(min_cells) || min_cells < 1) min_cells <- 30L
  writeLines(paste0("MIN_CELLS_PER_DONOR_SUBTYPE_CONDITION=", min_cells), project_path("outputs/tables/statistical_thresholds.txt"))

  neurons <- read_required_rds(project_path("outputs/objects/seurat_annotated_aplp1.rds"))
  neurons <- join_layers_if_needed(neurons, assay = "RNA")
  DefaultAssay(neurons) <- "RNA"

  counts <- get_counts_matrix(neurons, assay = "RNA")
  gene <- safe_feature("APLP1", counts)
  if (!length(gene)) stop("APLP1 not found in counts matrix.", call. = FALSE)

  donor_summary <- data.table::fread(project_path("outputs/tables/aplp1_fraction_by_donor_condition_subtype.csv"))
  donor_summary[, condition := as.character(condition)]
  donor_summary[, neuronal_subtype := as.character(neuronal_subtype)]

  subtypes <- sort(unique(donor_summary$neuronal_subtype))

  ci_ttest <- function(x, y = NULL, paired = FALSE) {
    out <- c(NA_real_, NA_real_)
    val <- tryCatch({
      if (paired) stats::t.test(x)$conf.int else stats::t.test(x, y)$conf.int
    }, error = function(e) NULL)
    if (!is.null(val)) out <- as.numeric(val[1:2])
    out
  }

  fraction_primary <- function(dt, subtype) {
    d <- dt[neuronal_subtype == subtype & condition %in% c("AD_AT8neg", "AD_AT8pos") & n_cells >= min_cells]
    pos <- d[condition == "AD_AT8pos"]
    neg <- d[condition == "AD_AT8neg"]
    wide <- merge(
      neg[, .(donor, AD_AT8neg = APLP1_positive_fraction)],
      pos[, .(donor, AD_AT8pos = APLP1_positive_fraction)],
      by = "donor"
    )
    if (nrow(wide) < 3) {
      return(empty_stats_row(subtype, "AD_AT8pos_vs_AD_AT8neg", "APLP1_positive_fraction", "insufficient paired donors or cells"))
    }
    delta <- wide$AD_AT8pos - wide$AD_AT8neg
    wt <- tryCatch(
      stats::wilcox.test(wide$AD_AT8pos, wide$AD_AT8neg, paired = TRUE, exact = FALSE),
      error = function(e) NULL
    )
    ci <- ci_ttest(delta, paired = TRUE)
    data.frame(
      neuronal_subtype = subtype,
      comparison = "AD_AT8pos_vs_AD_AT8neg",
      metric = "APLP1_positive_fraction",
      n_donors_group1 = data.table::uniqueN(pos$donor),
      n_donors_group2 = data.table::uniqueN(neg$donor),
      n_paired_donors = nrow(wide),
      effect_size = mean(delta),
      log2FC = NA_real_,
      ci_low = ci[1],
      ci_high = ci[2],
      p_value = if (!is.null(wt)) wt$p.value else NA_real_,
      FDR = NA_real_,
      status = "ok",
      stringsAsFactors = FALSE
    )
  }

  fraction_unpaired <- function(dt, subtype, group1, group2 = "nonAD") {
    comp <- paste0(group1, "_vs_", group2)
    d <- dt[neuronal_subtype == subtype & condition %in% c(group1, group2) & n_cells >= min_cells]
    x <- d[condition == group1]
    y <- d[condition == group2]
    if (data.table::uniqueN(x$donor) < 3 || data.table::uniqueN(y$donor) < 3) {
      return(empty_stats_row(subtype, comp, "APLP1_positive_fraction", "insufficient donors or cells"))
    }
    wt <- tryCatch(
      stats::wilcox.test(x$APLP1_positive_fraction, y$APLP1_positive_fraction, paired = FALSE, exact = FALSE),
      error = function(e) NULL
    )
    ci <- ci_ttest(x$APLP1_positive_fraction, y$APLP1_positive_fraction, paired = FALSE)
    data.frame(
      neuronal_subtype = subtype,
      comparison = comp,
      metric = "APLP1_positive_fraction",
      n_donors_group1 = data.table::uniqueN(x$donor),
      n_donors_group2 = data.table::uniqueN(y$donor),
      n_paired_donors = NA_integer_,
      effect_size = mean(x$APLP1_positive_fraction) - mean(y$APLP1_positive_fraction),
      log2FC = NA_real_,
      ci_low = ci[1],
      ci_high = ci[2],
      p_value = if (!is.null(wt)) wt$p.value else NA_real_,
      FDR = NA_real_,
      status = "ok",
      stringsAsFactors = FALSE
    )
  }

  fraction_stats <- data.table::rbindlist(c(
    lapply(subtypes, function(st) fraction_primary(donor_summary, st)),
    lapply(subtypes, function(st) fraction_unpaired(donor_summary, st, "AD_AT8pos")),
    lapply(subtypes, function(st) fraction_unpaired(donor_summary, st, "AD_AT8neg"))
  ), fill = TRUE)
  fraction_stats[
    !is.na(p_value),
    FDR := stats::p.adjust(p_value, method = "BH"),
    by = comparison
  ]
  write_table(fraction_stats, project_path("outputs/tables/aplp1_stats_fraction.csv"))

  md <- as.data.table(neurons@meta.data)
  md[, cell_barcode := rownames(neurons@meta.data)]
  md[, condition := as.character(condition)]
  md[, neuronal_subtype := as.character(neuronal_subtype)]
  md[, pb_group := paste(donor, condition, neuronal_subtype, sep = "|")]
  stopifnot(identical(colnames(counts), rownames(neurons@meta.data)))

  aggregate_subtype_counts <- function(subtype) {
    cells <- which(md$neuronal_subtype == subtype)
    if (!length(cells)) return(NULL)
    groups <- unique(md$pb_group[cells])
    pb <- do.call(cbind, lapply(groups, function(g) {
      idx <- cells[md$pb_group[cells] == g]
      Matrix::rowSums(counts[, idx, drop = FALSE])
    }))
    colnames(pb) <- groups
    rownames(pb) <- rownames(counts)
    group_meta <- data.table::rbindlist(lapply(groups, function(g) {
      parts <- strsplit(g, "\\|", fixed = FALSE)[[1]]
      idx <- cells[md$pb_group[cells] == g]
      data.table::data.table(
        sample = g,
        donor = parts[1],
        condition = parts[2],
        neuronal_subtype = parts[3],
        n_cells = length(idx)
      )
    }))
    list(counts = pb, meta = group_meta)
  }

  aplp1_pb_tables <- list()

  run_voom <- function(pb, pb_meta, subtype, comparison, group1, group2 = "nonAD", paired = FALSE) {
    d <- copy(pb_meta)
    d <- d[condition %in% c(group1, group2) & n_cells >= min_cells]
    if (paired) {
      wide_n <- dcast(d, donor ~ condition, value.var = "n_cells")
      if (!all(c(group1, group2) %in% colnames(wide_n))) {
        return(empty_stats_row(subtype, comparison, "APLP1_pseudobulk_expression", "insufficient paired donors or cells"))
      }
      wide_n <- wide_n[!is.na(get(group1)) & !is.na(get(group2))]
      donors <- wide_n$donor
      d <- d[donor %in% donors]
      if (length(donors) < 3) {
        return(empty_stats_row(subtype, comparison, "APLP1_pseudobulk_expression", "insufficient paired donors or cells"))
      }
    } else {
      if (data.table::uniqueN(d[condition == group1]$donor) < 3 || data.table::uniqueN(d[condition == group2]$donor) < 3) {
        return(empty_stats_row(subtype, comparison, "APLP1_pseudobulk_expression", "insufficient donors or cells"))
      }
    }
    sample_names <- d$sample
    y_counts <- pb[, sample_names, drop = FALSE]
    if (sum(y_counts[gene, , drop = TRUE]) == 0) {
      return(empty_stats_row(subtype, comparison, "APLP1_pseudobulk_expression", "APLP1 has zero pseudobulk counts in eligible samples"))
    }

    d[, condition := factor(condition, levels = c(group2, group1))]
    if (paired) {
      d[, donor := factor(donor)]
      design <- stats::model.matrix(~ donor + condition, data = d)
    } else {
      design <- stats::model.matrix(~ condition, data = d)
    }
    if (qr(design)$rank < ncol(design)) {
      return(empty_stats_row(subtype, comparison, "APLP1_pseudobulk_expression", "design matrix is not full rank"))
    }

    coef_name <- paste0("condition", group1)
    if (!coef_name %in% colnames(design)) {
      return(empty_stats_row(subtype, comparison, "APLP1_pseudobulk_expression", paste("missing coefficient", coef_name)))
    }

    y <- edgeR::DGEList(counts = y_counts)
    keep <- edgeR::filterByExpr(y, design = design)
    keep[gene] <- TRUE
    if (sum(keep) < 10) {
      return(empty_stats_row(subtype, comparison, "APLP1_pseudobulk_expression", "too few expressed genes for voom"))
    }
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    v <- limma::voom(y, design, plot = FALSE)
    fit <- limma::lmFit(v, design)
    fit <- limma::eBayes(fit)
    coef_idx <- which(colnames(design) == coef_name)
    if (!gene %in% rownames(fit$coefficients)) {
      return(empty_stats_row(subtype, comparison, "APLP1_pseudobulk_expression", "APLP1 absent after expression filtering"))
    }
    se <- fit$stdev.unscaled[gene, coef_idx] * sqrt(fit$s2.post[gene])
    lfc <- fit$coefficients[gene, coef_idx]
    p <- fit$p.value[gene, coef_idx]
    data.frame(
      neuronal_subtype = subtype,
      comparison = comparison,
      metric = "APLP1_pseudobulk_expression",
      n_donors_group1 = data.table::uniqueN(d[condition == group1]$donor),
      n_donors_group2 = data.table::uniqueN(d[condition == group2]$donor),
      n_paired_donors = if (paired) data.table::uniqueN(d$donor) else NA_integer_,
      effect_size = NA_real_,
      log2FC = lfc,
      ci_low = lfc - 1.96 * se,
      ci_high = lfc + 1.96 * se,
      p_value = p,
      FDR = NA_real_,
      status = "ok",
      stringsAsFactors = FALSE
    )
  }

  pb_stats_list <- list()
  for (st in subtypes) {
    agg <- aggregate_subtype_counts(st)
    if (is.null(agg)) next
    pb <- agg$counts
    pb_meta <- agg$meta
    lib_size <- Matrix::colSums(pb)
    aplp1_count <- as.numeric(pb[gene, , drop = TRUE])
    aplp1_pb_tables[[st]] <- cbind(
      pb_meta,
      data.table::data.table(
        library_size = as.numeric(lib_size),
        APLP1_pseudobulk_count = aplp1_count,
        APLP1_log2CPM = log2(((aplp1_count + 0.5) / (as.numeric(lib_size) + 1)) * 1e6)
      )
    )

    pb_stats_list[[paste(st, "primary")]] <- run_voom(
      pb, pb_meta, st, "AD_AT8pos_vs_AD_AT8neg", "AD_AT8pos", "AD_AT8neg", paired = TRUE
    )
    pb_stats_list[[paste(st, "pos_nonad")]] <- run_voom(
      pb, pb_meta, st, "AD_AT8pos_vs_nonAD", "AD_AT8pos", "nonAD", paired = FALSE
    )
    pb_stats_list[[paste(st, "neg_nonad")]] <- run_voom(
      pb, pb_meta, st, "AD_AT8neg_vs_nonAD", "AD_AT8neg", "nonAD", paired = FALSE
    )
  }

  aplp1_pb_table <- data.table::rbindlist(aplp1_pb_tables, fill = TRUE)
  write_table(aplp1_pb_table, project_path("outputs/tables/aplp1_pseudobulk_by_donor_condition_subtype.csv"))

  pb_stats <- data.table::rbindlist(pb_stats_list, fill = TRUE)
  pb_stats[
    !is.na(p_value),
    FDR := stats::p.adjust(p_value, method = "BH"),
    by = comparison
  ]
  write_table(pb_stats, project_path("outputs/tables/aplp1_stats_pseudobulk.csv"))

  all_stats <- data.table::rbindlist(list(fraction_stats, pb_stats), fill = TRUE)
  write_table(all_stats, project_path("outputs/tables/aplp1_stats_all.csv"))

  primary_summary <- all_stats[comparison == "AD_AT8pos_vs_AD_AT8neg"][
    order(metric, FDR, p_value, neuronal_subtype)
  ]
  write_table(primary_summary, project_path("outputs/tables/aplp1_primary_comparison_summary.csv"))
  log_message("Wrote donor-level fraction and pseudobulk statistics with min_cells =", min_cells, log = "05_aplp1_stats.log")
})

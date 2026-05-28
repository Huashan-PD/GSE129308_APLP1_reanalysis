#!/usr/bin/env Rscript

source("scripts/lib_project.R")

run_main("03_cluster_annotate", {
  load_package("Seurat", required = TRUE, log = "03_cluster_annotate.log")
  load_package("SeuratObject", required = TRUE, log = "03_cluster_annotate.log")
  load_package("Matrix", required = TRUE, log = "03_cluster_annotate.log")
  load_package("data.table", required = TRUE, log = "03_cluster_annotate.log")
  load_package("ggplot2", required = TRUE, log = "03_cluster_annotate.log")
  load_package("patchwork", required = FALSE, log = "03_cluster_annotate.log")
  load_package("scales", required = FALSE, log = "03_cluster_annotate.log")

  options(future.globals.maxSize = 50 * 1024^3)

  obj <- read_required_rds(project_path("outputs/objects/seurat_qc.rds"))
  obj <- join_layers_if_needed(obj, assay = "RNA")
  DefaultAssay(obj) <- "RNA"
  obj$condition <- as_condition_factor(obj$condition)

  run_workflow <- function(x, label, resolution = 0.6, nfeatures = 3000, dims = 1:30) {
    set.seed(PROJECT_SEED)
    x <- join_layers_if_needed(x, assay = "RNA")
    DefaultAssay(x) <- "RNA"
    x <- Seurat::NormalizeData(x, verbose = FALSE)
    x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)

    use_cca <- tolower(Sys.getenv("USE_CCA_INTEGRATION", unset = "1")) %in% c("1", "true", "yes")
    method_used <- "LogNormalize plus PCA without sample integration"
    reduction_use <- "pca"
    integrated_done <- FALSE

    integration_method <- tolower(Sys.getenv("INTEGRATION_METHOD", unset = "harmony"))
    if (use_cca && integration_method == "harmony" && requireNamespace("harmony", quietly = TRUE)) {
      harmony_attempt <- tryCatch(
        {
          log_message("Attempting Harmony integration for", label, "by sample_id.", log = "03_cluster_annotate.log")
          vars_to_regress <- intersect(c("nCount_RNA", "percent.mt"), colnames(x@meta.data))
          x <- Seurat::ScaleData(x, vars.to.regress = vars_to_regress, verbose = FALSE)
          x <- Seurat::RunPCA(x, npcs = max(dims), verbose = FALSE)
          x <- harmony::RunHarmony(
            x,
            group.by.vars = "sample_id",
            reduction.use = "pca",
            dims.use = dims,
            reduction.save = "harmony",
            verbose = FALSE
          )
          x
        },
        error = function(e) {
          log_message("Harmony integration failed for", label, ":", conditionMessage(e), log = "03_cluster_annotate.log")
          NULL
        }
      )
      if (!is.null(harmony_attempt)) {
        x <- harmony_attempt
        method_used <- "Harmony integration by sample_id"
        reduction_use <- "harmony"
        integrated_done <- TRUE
      }
    }

    integration_reduction <- tolower(Sys.getenv("SEURAT_INTEGRATION_REDUCTION", unset = "rpca"))
    if (!integration_reduction %in% c("rpca", "cca")) integration_reduction <- "rpca"

    if (!integrated_done && use_cca && length(unique(x$sample_id)) > 1) {
      integrated_attempt <- tryCatch(
        {
          log_message(
            "Attempting Seurat", toupper(integration_reduction), "integration for", label, "by sample_id.",
            log = "03_cluster_annotate.log"
          )
          split_obj <- Seurat::SplitObject(x, split.by = "sample_id")
          split_obj <- lapply(split_obj, function(y) {
            y <- Seurat::NormalizeData(y, verbose = FALSE)
            y <- Seurat::FindVariableFeatures(y, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
            y
          })
          features <- Seurat::SelectIntegrationFeatures(object.list = split_obj, nfeatures = nfeatures)
          if (integration_reduction == "rpca") {
            split_obj <- lapply(split_obj, function(y) {
              y <- Seurat::ScaleData(y, features = features, verbose = FALSE)
              y <- Seurat::RunPCA(y, features = features, npcs = max(dims), verbose = FALSE)
              y
            })
          }
          anchors <- Seurat::FindIntegrationAnchors(
            object.list = split_obj,
            anchor.features = features,
            dims = dims,
            reduction = integration_reduction,
            verbose = FALSE
          )
          z <- Seurat::IntegrateData(anchorset = anchors, dims = dims, verbose = FALSE)
          DefaultAssay(z) <- "integrated"
          z <- Seurat::ScaleData(z, verbose = FALSE)
          z <- Seurat::RunPCA(z, npcs = max(dims), verbose = FALSE)
          method_used <- paste0("Seurat ", toupper(integration_reduction), " integration by sample_id")
          z
        },
        error = function(e) {
          log_message("CCA integration failed for", label, ":", conditionMessage(e), log = "03_cluster_annotate.log")
          NULL
        }
      )
      if (!is.null(integrated_attempt)) x <- integrated_attempt
    }

    if (!"pca" %in% names(x@reductions)) {
      DefaultAssay(x) <- "RNA"
      vars_to_regress <- intersect(c("nCount_RNA", "percent.mt"), colnames(x@meta.data))
      x <- Seurat::ScaleData(x, vars.to.regress = vars_to_regress, verbose = FALSE)
      x <- Seurat::RunPCA(x, npcs = max(dims), verbose = FALSE)
    }

    if (!reduction_use %in% names(x@reductions)) reduction_use <- "pca"
    ndims <- min(max(dims), ncol(Seurat::Embeddings(x, reduction_use)))
    dims_use <- seq_len(ndims)
    x <- Seurat::RunUMAP(x, reduction = reduction_use, dims = dims_use, seed.use = PROJECT_SEED, verbose = FALSE)
    x <- Seurat::FindNeighbors(x, reduction = reduction_use, dims = dims_use, verbose = FALSE)
    x <- Seurat::FindClusters(x, resolution = resolution, random.seed = PROJECT_SEED, verbose = FALSE)
    x$integration_method <- method_used
    x
  }

  score_marker_sets <- function(x, marker_sets, group_by, min_score = 0.05) {
    DefaultAssay(x) <- "RNA"
    features <- unique(unlist(marker_sets, use.names = FALSE))
    features <- present_features(features, x)
    if (!length(features)) stop("None of the requested marker genes are present in the object.", call. = FALSE)
    avg <- Seurat::AverageExpression(
      x,
      assays = "RNA",
      features = features,
      group.by = group_by,
      slot = "data",
      verbose = FALSE
    )$RNA
    groups_raw <- colnames(avg)
    groups <- sub("^g(?=[0-9])", "", groups_raw, perl = TRUE)
    score_list <- lapply(names(marker_sets), function(label) {
      markers <- intersect(marker_sets[[label]], rownames(avg))
      if (!length(markers)) {
        score <- rep(NA_real_, length(groups))
      } else {
        score <- Matrix::colMeans(avg[markers, , drop = FALSE])
      }
      data.frame(group = groups, annotation = label, score = score, markers_present = paste(markers, collapse = ";"), stringsAsFactors = FALSE)
    })
    scores <- data.table::rbindlist(score_list)
    assignments <- scores[
      ,
      {
        valid <- .SD[!is.na(score)]
        if (!nrow(valid) || max(valid$score) < min_score) {
          data.table::data.table(annotation = "Unassigned", score = ifelse(nrow(valid), max(valid$score), NA_real_))
        } else {
          best <- valid[which.max(score)]
          data.table::data.table(annotation = best$annotation, score = best$score)
        }
      },
      by = group
    ]
    list(scores = as.data.frame(scores), assignments = as.data.frame(assignments))
  }

  cached_all_cells <- project_path("outputs/objects/seurat_annotated_all_cells.rds")
  force_recluster_all <- tolower(Sys.getenv("FORCE_RECLUSTER_ALL", unset = "0")) %in% c("1", "true", "yes")
  if (file.exists(cached_all_cells) && !force_recluster_all) {
    log_message("Reusing cached all-cell clustering object:", cached_all_cells, log = "03_cluster_annotate.log")
    obj <- readRDS(cached_all_cells)
    obj <- join_layers_if_needed(obj, assay = "RNA")
    DefaultAssay(obj) <- "RNA"
    if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
      stop("Cached all-cell object lacks seurat_clusters; set FORCE_RECLUSTER_ALL=1.", call. = FALSE)
    }
  } else {
    obj <- run_workflow(obj, label = "all QC cells", resolution = 0.6)
  }

  major_markers <- list(
    Excitatory_neuron = c("SLC17A7"),
    Inhibitory_neuron = c("GAD1", "GAD2"),
    Astrocyte = c("AQP4", "GFAP", "SLC1A3"),
    Oligodendrocyte = c("MBP", "MOG", "PLP1"),
    OPC = c("PDGFRA", "CSPG4"),
    Microglia = c("CX3CR1", "P2RY12", "C1QA"),
    Endothelial = c("CLDN5", "FLT1")
  )
  major_scores <- score_marker_sets(obj, major_markers, group_by = "seurat_clusters", min_score = 0.05)
  write_table(major_scores$scores, project_path("outputs/tables/major_celltype_cluster_marker_scores.csv"))
  write_table(major_scores$assignments, project_path("outputs/tables/major_celltype_cluster_assignments.csv"))

  cluster_to_major <- setNames(major_scores$assignments$annotation, major_scores$assignments$group)
  obj$major_cell_type <- unname(cluster_to_major[as.character(obj$seurat_clusters)])
  obj$major_cell_type[is.na(obj$major_cell_type)] <- "Unassigned"

  save_rds(obj, project_path("outputs/objects/seurat_annotated_all_cells.rds"))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    save_plot_both(Seurat::DimPlot(obj, group.by = "condition", reduction = "umap"), "umap_allcells_condition", width = 7, height = 5)
    save_plot_both(Seurat::DimPlot(obj, group.by = "donor", reduction = "umap"), "umap_allcells_donor", width = 8, height = 5)
    save_plot_both(Seurat::DimPlot(obj, group.by = "sample_id", reduction = "umap"), "umap_allcells_sample", width = 10, height = 6)
    save_plot_both(Seurat::DimPlot(obj, group.by = "major_cell_type", reduction = "umap", label = TRUE, repel = TRUE), "umap_allcells_major_cell_type", width = 8, height = 6)
  }

  neuronal_cells <- colnames(obj)[obj$major_cell_type %in% c("Excitatory_neuron", "Inhibitory_neuron")]
  if (length(neuronal_cells) < 100) {
    stop("Fewer than 100 neuronal cells were annotated. Check marker mapping and major cell type assignments.", call. = FALSE)
  }

  neurons <- subset(obj, cells = neuronal_cells)
  DefaultAssay(neurons) <- "RNA"
  neurons <- run_workflow(neurons, label = "neuronal cells", resolution = 0.8)
  neurons$neuronal_cluster <- as.character(neurons$seurat_clusters)

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

  subtype_scores <- score_marker_sets(neurons, subtype_markers, group_by = "neuronal_cluster", min_score = 0.03)
  write_table(subtype_scores$scores, project_path("outputs/tables/neuronal_subtype_cluster_marker_scores.csv"))

  cluster_major_mode <- as.data.table(neurons@meta.data)[
    ,
    .N,
    by = .(neuronal_cluster, major_cell_type)
  ][
    order(neuronal_cluster, -N)
  ][
    ,
    .SD[1],
    by = neuronal_cluster
  ]

  assignments <- as.data.table(subtype_scores$assignments)
  assignments <- merge(assignments, cluster_major_mode, by.x = "group", by.y = "neuronal_cluster", all.x = TRUE)
  assignments[
    major_cell_type == "Excitatory_neuron" & !grepl("^Ex", annotation),
    annotation := "Excitatory_unresolved"
  ]
  assignments[
    major_cell_type == "Inhibitory_neuron" & !grepl("^In", annotation),
    annotation := "Inhibitory_unresolved"
  ]
  setnames(assignments, "group", "neuronal_cluster")
  write_table(assignments, project_path("outputs/tables/neuronal_subtype_cluster_assignments.csv"))

  cluster_to_subtype <- setNames(assignments$annotation, assignments$neuronal_cluster)
  neurons$neuronal_subtype <- unname(cluster_to_subtype[as.character(neurons$neuronal_cluster)])
  neurons$neuronal_subtype[is.na(neurons$neuronal_subtype)] <- "Unassigned"
  neurons$condition <- as_condition_factor(neurons$condition)

  cell_meta <- as.data.frame(neurons@meta.data)
  cell_meta$cell_barcode <- rownames(cell_meta)
  write_table_gz(cell_meta, project_path("outputs/tables/cell_metadata_annotated.csv.gz"))

  subtype_counts <- as.data.table(neurons@meta.data)[
    ,
    .(n_cells = .N, n_donors = data.table::uniqueN(donor), n_samples = data.table::uniqueN(sample_id)),
    by = .(neuronal_subtype, condition)
  ][order(neuronal_subtype, condition)]
  write_table(subtype_counts, project_path("outputs/tables/neuronal_subtype_cell_counts.csv"))

  marker_features <- present_features(unique(unlist(subtype_markers, use.names = FALSE)), neurons)
  if (length(marker_features)) {
    p_marker <- Seurat::DotPlot(
      neurons,
      features = marker_features,
      group.by = "neuronal_subtype",
      assay = "RNA"
    ) + ggplot2::coord_flip() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7))
    save_plot_both(p_marker, "marker_dotplot_neuronal_subtypes", width = 11, height = 8)
  }

  save_plot_both(Seurat::DimPlot(neurons, group.by = "condition", reduction = "umap"), "umap_neurons_condition", width = 7, height = 5)
  save_plot_both(Seurat::DimPlot(neurons, group.by = "donor", reduction = "umap"), "umap_neurons_donor", width = 8, height = 5)
  save_plot_both(Seurat::DimPlot(neurons, group.by = "sample_id", reduction = "umap"), "umap_neurons_sample", width = 10, height = 6)
  save_plot_both(
    Seurat::DimPlot(neurons, group.by = "neuronal_subtype", reduction = "umap", label = TRUE, repel = TRUE),
    "umap_neurons_neuronal_subtype",
    width = 10,
    height = 7
  )

  writeLines(unique(neurons$integration_method), project_path("outputs/tables/integration_method.txt"))
  save_rds(neurons, project_path("outputs/objects/seurat_annotated.rds"))
  log_message("Annotated neuronal object contains", ncol(neurons), "cells.", log = "03_cluster_annotate.log")
})

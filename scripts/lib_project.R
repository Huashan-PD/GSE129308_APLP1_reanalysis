PROJECT_SEED <- 20260526
set.seed(PROJECT_SEED)

find_project_root <- function() {
  env_root <- Sys.getenv("PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = TRUE))
  }
  wd <- normalizePath(getwd(), mustWork = TRUE)
  if (basename(wd) == "scripts" && file.exists(file.path(wd, "lib_project.R"))) {
    return(dirname(wd))
  }
  if (basename(wd) == "GSE129308_APLP1_reanalysis" && dir.exists(file.path(wd, "scripts"))) {
    return(wd)
  }
  candidate <- file.path(wd, "GSE129308_APLP1_reanalysis")
  if (dir.exists(file.path(candidate, "scripts"))) {
    return(normalizePath(candidate, mustWork = TRUE))
  }
  stop("Could not locate project root. Run scripts from GSE129308_APLP1_reanalysis or set PROJECT_ROOT.")
}

PROJECT_ROOT <- find_project_root()

project_path <- function(...) file.path(PROJECT_ROOT, ...)

local_r_lib <- project_path("r_libs")
dir.create(local_r_lib, recursive = TRUE, showWarnings = FALSE)
if (dir.exists(local_r_lib)) {
  .libPaths(unique(c(local_r_lib, .libPaths())))
}

project_dirs <- c(
  "data/raw", "data/intermediate", "data/processed",
  "scripts", "outputs/figures", "outputs/tables",
  "outputs/objects", "outputs/report", "logs"
)

ensure_project_dirs <- function() {
  for (d in project_dirs) {
    dir.create(project_path(d), recursive = TRUE, showWarnings = FALSE)
  }
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

log_message <- function(..., log = "pipeline.log") {
  ensure_project_dirs()
  msg <- paste0("[", timestamp(), "] ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = project_path("logs", log), append = TRUE)
}

log_error <- function(..., script = NULL) {
  ensure_project_dirs()
  prefix <- if (!is.null(script)) paste0("[", script, "] ") else ""
  msg <- paste0("[", timestamp(), "] ERROR ", prefix, paste(..., collapse = " "))
  cat(msg, "\n", file = stderr())
  cat(msg, "\n", file = project_path("logs", "error.log"), append = TRUE)
}

run_main <- function(script_name, expr) {
  ensure_project_dirs()
  log_message("Starting", script_name, log = paste0(script_name, ".log"))
  tryCatch(
    {
      force(expr)
      log_message("Finished", script_name, log = paste0(script_name, ".log"))
    },
    error = function(e) {
      log_error(conditionMessage(e), script = script_name)
      log_message("Failed", script_name, conditionMessage(e), log = paste0(script_name, ".log"))
      quit(status = 1, save = "no")
    }
  )
}

load_package <- function(pkg, required = TRUE, log = "pipeline.log") {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok) {
    msg <- paste0("R package not available: ", pkg)
    log_message(msg, log = log)
    if (required) stop(msg, call. = FALSE)
    return(FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  TRUE
}

write_table <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, path, ...)
  } else {
    utils::write.csv(x, path, row.names = FALSE, ...)
  }
}

write_table_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE)
}

save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path)
}

read_required_rds <- function(path) {
  if (!file.exists(path)) stop("Required object is missing: ", path, call. = FALSE)
  readRDS(path)
}

save_plot_both <- function(plot, filename_base, width = 8, height = 6, dpi = 300) {
  load_package("ggplot2", required = TRUE)
  dir.create(project_path("outputs/figures"), recursive = TRUE, showWarnings = FALSE)
  pdf_path <- project_path("outputs/figures", paste0(filename_base, ".pdf"))
  png_path <- project_path("outputs/figures", paste0(filename_base, ".png"))
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
  invisible(c(pdf = pdf_path, png = png_path))
}

condition_levels <- c("nonAD", "AD_AT8neg", "AD_AT8pos")
condition_labels <- c(
  nonAD = "non-AD control neurons/somas",
  AD_AT8neg = "AD-AT8- neighboring neurons/somas",
  AD_AT8pos = "AD-AT8+ NFT-bearing neurons/somas"
)

condition_short_labels <- c(
  nonAD = "nonAD",
  AD_AT8neg = "AD_AT8neg",
  AD_AT8pos = "AD_AT8pos"
)

as_condition_factor <- function(x) factor(x, levels = condition_levels)

get_counts_matrix <- function(obj, assay = "RNA") {
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    out <- tryCatch(
      SeuratObject::LayerData(object = obj, assay = assay, layer = "counts"),
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
    return(SeuratObject::GetAssayData(object = obj, assay = assay, slot = "counts"))
  }
  stop("SeuratObject is required to extract counts.")
}

get_data_matrix <- function(obj, assay = "RNA") {
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    out <- tryCatch(
      SeuratObject::LayerData(object = obj, assay = assay, layer = "data"),
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
    return(SeuratObject::GetAssayData(object = obj, assay = assay, slot = "data"))
  }
  stop("SeuratObject is required to extract normalized data.")
}

join_layers_if_needed <- function(obj, assay = "RNA") {
  if (requireNamespace("SeuratObject", quietly = TRUE) &&
      exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
    obj <- tryCatch(
      SeuratObject::JoinLayers(obj, assay = assay),
      error = function(e) obj
    )
  }
  obj
}

present_features <- function(features, object_or_matrix) {
  rn <- rownames(object_or_matrix)
  features[features %in% rn]
}

safe_feature <- function(feature, object_or_matrix) {
  rn <- rownames(object_or_matrix)
  if (feature %in% rn) return(feature)
  hit <- rn[toupper(rn) == toupper(feature)]
  if (length(hit) == 1) return(hit)
  character(0)
}

empty_stats_row <- function(subtype, comparison, metric, reason) {
  data.frame(
    neuronal_subtype = subtype,
    comparison = comparison,
    metric = metric,
    n_donors_group1 = NA_integer_,
    n_donors_group2 = NA_integer_,
    n_paired_donors = NA_integer_,
    effect_size = NA_real_,
    log2FC = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    p_value = NA_real_,
    FDR = NA_real_,
    status = reason,
    stringsAsFactors = FALSE
  )
}

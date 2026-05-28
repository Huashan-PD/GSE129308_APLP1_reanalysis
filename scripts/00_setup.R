#!/usr/bin/env Rscript

source("scripts/lib_project.R")

run_main("00_setup", {
  required_pkgs <- c(
    "Seurat", "SeuratObject", "Matrix", "hdf5r", "tidyverse", "data.table",
    "ggplot2", "patchwork", "pheatmap", "ComplexHeatmap", "edgeR", "limma",
    "glmmTMB", "betareg", "lme4", "MAST", "cowplot", "ggrepel", "scales",
    "knitr", "rmarkdown"
  )
  bioc_pkgs <- c("ComplexHeatmap", "edgeR", "limma", "MAST")
  cran_pkgs <- setdiff(required_pkgs, bioc_pkgs)

  install_missing <- tolower(Sys.getenv("INSTALL_MISSING_R", unset = "0")) %in% c("1", "true", "yes")
  installed_before <- rownames(installed.packages())
  missing_before <- setdiff(required_pkgs, installed_before)

  if (length(missing_before)) {
    log_message("Missing packages before setup:", paste(missing_before, collapse = ", "), log = "00_setup.log")
  } else {
    log_message("All requested packages are already installed.", log = "00_setup.log")
  }

  if (install_missing && length(missing_before)) {
    for (pkg in intersect(missing_before, cran_pkgs)) {
      log_message("Attempting CRAN install:", pkg, log = "00_setup.log")
      tryCatch(
        utils::install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE),
        error = function(e) log_message("CRAN install failed for", pkg, ":", conditionMessage(e), log = "00_setup.log")
      )
    }
    if (length(intersect(missing_before, bioc_pkgs))) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        tryCatch(
          utils::install.packages("BiocManager", repos = "https://cloud.r-project.org"),
          error = function(e) log_message("BiocManager install failed:", conditionMessage(e), log = "00_setup.log")
        )
      }
      if (requireNamespace("BiocManager", quietly = TRUE)) {
        for (pkg in intersect(missing_before, bioc_pkgs)) {
          log_message("Attempting Bioconductor install:", pkg, log = "00_setup.log")
          tryCatch(
            BiocManager::install(pkg, ask = FALSE, update = FALSE),
            error = function(e) log_message("Bioconductor install failed for", pkg, ":", conditionMessage(e), log = "00_setup.log")
          )
        }
      }
    }
  } else if (length(missing_before)) {
    log_message(
      "Set INSTALL_MISSING_R=1 before running 00_setup.R to attempt package installation.",
      log = "00_setup.log"
    )
  }

  installed_after <- rownames(installed.packages())
  status <- data.frame(
    package = required_pkgs,
    installed = required_pkgs %in% installed_after,
    required_for_core = required_pkgs %in% c(
      "Seurat", "SeuratObject", "Matrix", "hdf5r", "data.table", "ggplot2",
      "patchwork", "edgeR", "limma", "knitr", "rmarkdown"
    ),
    stringsAsFactors = FALSE
  )
  write_table(status, project_path("outputs/tables/package_status.csv"))
  missing_after <- status$package[!status$installed]
  writeLines(missing_after, project_path("outputs/tables/missing_packages.txt"))

  sink(project_path("outputs/tables/sessionInfo.txt"))
  print(sessionInfo())
  sink()

  if (length(missing_after)) {
    log_message("Packages still missing:", paste(missing_after, collapse = ", "), log = "00_setup.log")
  }
})


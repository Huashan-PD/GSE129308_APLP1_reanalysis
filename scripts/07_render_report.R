#!/usr/bin/env Rscript

source("scripts/lib_project.R")

run_main("07_render_report", {
  load_package("rmarkdown", required = TRUE, log = "07_render_report.log")
  load_package("knitr", required = TRUE, log = "07_render_report.log")

  rmd <- project_path("outputs/report/APLP1_GSE129308_report.Rmd")
  if (!file.exists(rmd)) stop("Report Rmd is missing: ", rmd, call. = FALSE)

  if (!rmarkdown::pandoc_available()) {
    pandoc_candidates <- c(
      "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/pandoc",
      "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64/pandoc",
      "/Applications/RStudio.app/Contents/MacOS/pandoc/pandoc"
    )
    pandoc_hit <- pandoc_candidates[file.exists(pandoc_candidates)][1]
    if (!is.na(pandoc_hit)) {
      Sys.setenv(RSTUDIO_PANDOC = dirname(pandoc_hit))
      log_message("Using pandoc at:", pandoc_hit, log = "07_render_report.log")
    }
  }
  if (!rmarkdown::pandoc_available()) {
    stop("pandoc is required for report rendering and was not found.", call. = FALSE)
  }

  html_out <- rmarkdown::render(
    input = rmd,
    output_format = "html_document",
    output_file = "APLP1_GSE129308_report.html",
    output_dir = project_path("outputs/report"),
    knit_root_dir = PROJECT_ROOT,
    envir = new.env(parent = globalenv()),
    quiet = FALSE
  )
  log_message("Rendered HTML report:", html_out, log = "07_render_report.log")

  pdf_out <- tryCatch(
    rmarkdown::render(
      input = rmd,
      output_format = "pdf_document",
      output_file = "APLP1_GSE129308_report.pdf",
      output_dir = project_path("outputs/report"),
      knit_root_dir = PROJECT_ROOT,
      envir = new.env(parent = globalenv()),
      quiet = FALSE
    ),
    error = function(e) {
      log_message("PDF report was not rendered:", conditionMessage(e), log = "07_render_report.log")
      NULL
    }
  )
  if (!is.null(pdf_out)) log_message("Rendered PDF report:", pdf_out, log = "07_render_report.log")
})

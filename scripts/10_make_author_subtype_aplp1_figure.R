#!/usr/bin/env Rscript

source("scripts/lib_project.R")

run_main("10_make_author_subtype_aplp1_figure", {
  load_package("data.table", required = TRUE, log = "10_make_author_subtype_aplp1_figure.log")
  load_package("ggplot2", required = TRUE, log = "10_make_author_subtype_aplp1_figure.log")
  load_package("scales", required = TRUE, log = "10_make_author_subtype_aplp1_figure.log")
  load_package("patchwork", required = TRUE, log = "10_make_author_subtype_aplp1_figure.log")

  table_s2 <- project_path("data/raw/author_supplement/NIHMS1821336-supplement-Table_S2.xlsx")
  table_s6 <- project_path("data/raw/author_supplement/NIHMS1821336-supplement-Table_S6.xlsx")
  if (!file.exists(table_s2)) {
    stop("Author supplementary Table S2 is missing: ", table_s2, call. = FALSE)
  }
  if (!file.exists(table_s6)) {
    stop("Author supplementary Table S6 is missing: ", table_s6, call. = FALSE)
  }

  have_readxl <- requireNamespace("readxl", quietly = TRUE)
  cached_aplp1 <- project_path("outputs/tables/author_tableS6_APLP1_rows.csv")
  if (!have_readxl && !file.exists(cached_aplp1)) {
    stop("Package readxl is required to parse author xlsx files, and no cached APLP1 table was found.", call. = FALSE)
  }

  metadata_dir <- project_path("outputs/tables/author_metadata_search")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

  numeric_clean <- function(x) {
    suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
  }

  if (have_readxl) {
    inventory <- data.table::rbindlist(lapply(c(table_s2, table_s6), function(path) {
      sheets <- readxl::excel_sheets(path)
      data.table::rbindlist(lapply(sheets, function(sheet) {
        preview <- suppressWarnings(readxl::read_excel(path, sheet = sheet, n_max = 50))
        names_preview <- names(preview)
        barcode_like <- any(grepl(
          "(^|_|\\b)(barcode|cell[_ ]?barcode|cell[_ ]?id|soma[_ ]?id)(_|\\b|$)",
          names_preview,
          ignore.case = TRUE
        ))
        data.table::data.table(
          workbook = basename(path),
          sheet = sheet,
          preview_rows = nrow(preview),
          preview_cols = ncol(preview),
          barcode_like_column_detected = barcode_like,
          columns = paste(names_preview, collapse = " | ")
        )
      }), fill = TRUE)
    }), fill = TRUE)
    write_table(inventory, file.path(metadata_dir, "author_supplement_inventory.csv"))

    summary_s6 <- suppressWarnings(readxl::read_excel(table_s6, sheet = "SUMMARY", col_names = FALSE))
    summary_s6 <- data.table::as.data.table(summary_s6)
    data.table::setnames(summary_s6, paste0("V", seq_len(ncol(summary_s6))))
    subtype_counts <- summary_s6[grepl("^(Ex|In)[0-9]+\\s*\\(", V1)]
    subtype_counts[, subtype := sub("^((Ex|In)[0-9]+).*", "\\1", V1)]
    subtype_counts <- subtype_counts[
      ,
      .(
        subtype,
        subtype_description = V1,
        n_total = as.integer(numeric_clean(V2)),
        n_nonAD = as.integer(numeric_clean(V3)),
        n_AD_AT8neg = as.integer(numeric_clean(V4)),
        n_AD_AT8pos = as.integer(numeric_clean(V5))
      )
    ]
    write_table(subtype_counts, project_path("outputs/tables/author_tableS6_subtype_counts.csv"))

    s6_sheets <- setdiff(readxl::excel_sheets(table_s6), "SUMMARY")
    aplp1_rows <- data.table::rbindlist(lapply(s6_sheets, function(sheet) {
      dt <- suppressWarnings(data.table::as.data.table(readxl::read_excel(table_s6, sheet = sheet)))
      gene_col <- names(dt)[tolower(names(dt)) == "gene"][1]
      if (is.na(gene_col)) {
        return(NULL)
      }
      hit <- dt[toupper(as.character(get(gene_col))) == "APLP1"]
      if (!nrow(hit)) {
        return(NULL)
      }
      hit[, subtype := sheet]
      hit
    }), fill = TRUE)
    if (!nrow(aplp1_rows)) {
      stop("APLP1 was not found in author supplementary Table S6.", call. = FALSE)
    }
    write_table(aplp1_rows, cached_aplp1)
  } else {
    aplp1_rows <- data.table::fread(cached_aplp1)
    subtype_counts <- data.table::fread(project_path("outputs/tables/author_tableS6_subtype_counts.csv"))
  }

  required_cols <- c("pct.1 (AD-AT8+)", "pct.2 (non-AD)", "pct.3 (AD-AT8-)")
  missing_cols <- setdiff(required_cols, names(aplp1_rows))
  if (length(missing_cols)) {
    stop("Required Table S6 APLP1 pct columns are missing: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  primary_p_col <- "AD-AT8+ vs AD-AT8- p_val_adj"

  row_order <- c("Ex1", "Ex2", "Ex3", "Ex4&5", "Ex6", "Ex7", "Ex8", "Ex10", "Ex13", "In2")
  conditions <- data.table::data.table(
    condition_plot = c("AD-AT8+", "non-AD control"),
    pct_col = c("pct.1 (AD-AT8+)", "pct.2 (non-AD)"),
    count_col = c("n_AD_AT8pos", "n_nonAD")
  )

  get_weighted_fraction <- function(source_subtypes, pct_col, count_col) {
    x <- aplp1_rows[subtype %in% source_subtypes, .(subtype, fraction = numeric_clean(get(pct_col)))]
    w <- subtype_counts[subtype %in% source_subtypes, .(subtype, n_cells = numeric_clean(get(count_col)))]
    merged <- merge(x, w, by = "subtype", all.x = TRUE, sort = FALSE)
    if (!nrow(merged)) {
      return(list(fraction = NA_real_, n_cells = NA_real_, method = "missing Table S6 value"))
    }
    n_total <- if (all(!is.na(merged$n_cells)) && sum(merged$n_cells) > 0) {
      sum(merged$n_cells)
    } else {
      NA_real_
    }
    if (anyNA(merged$fraction)) {
      return(list(fraction = NA_real_, n_cells = n_total, method = "missing Table S6 value"))
    }
    if (!is.na(n_total)) {
      return(list(
        fraction = stats::weighted.mean(merged$fraction, merged$n_cells),
        n_cells = n_total,
        method = if (length(source_subtypes) > 1) {
          "cell-count weighted mean of author Table S6 subtypes"
        } else {
          "direct author Table S6 pct value"
        }
      ))
    }
    list(
      fraction = mean(merged$fraction),
      n_cells = NA_real_,
      method = if (length(source_subtypes) > 1) {
        "unweighted mean of author Table S6 subtypes; counts unavailable"
      } else {
        "direct author Table S6 pct value; counts unavailable"
      }
    )
  }

  fig_dt <- data.table::rbindlist(lapply(row_order, function(display_subtype) {
    source_subtypes <- if (display_subtype == "Ex4&5") c("Ex4", "Ex5") else display_subtype
    data.table::rbindlist(lapply(seq_len(nrow(conditions)), function(i) {
      res <- get_weighted_fraction(source_subtypes, conditions$pct_col[i], conditions$count_col[i])
      data.table::data.table(
        neuronal_subtype = display_subtype,
        source_subtypes = paste(source_subtypes, collapse = "+"),
        condition = conditions$condition_plot[i],
        APLP1_positive_fraction = res$fraction,
        author_n_cells = res$n_cells,
        source_column = conditions$pct_col[i],
        combination_method = res$method
      )
    }))
  }))

  fig_dt[, neuronal_subtype := factor(neuronal_subtype, levels = rev(row_order))]
  fig_dt[, condition := factor(condition, levels = c("AD-AT8+", "non-AD control"))]
  write_table(fig_dt, project_path("outputs/tables/author_tableS6_APLP1_fraction_author_subtypes.csv"))

  wide_dt <- data.table::dcast(
    fig_dt,
    neuronal_subtype + source_subtypes ~ condition,
    value.var = "APLP1_positive_fraction"
  )
  wide_dt[, neuronal_subtype := as.character(neuronal_subtype)]
  wide_dt[, row_order_i := match(neuronal_subtype, row_order)]
  data.table::setorder(wide_dt, row_order_i)
  wide_dt[, row_order_i := NULL]
  write_table(wide_dt, project_path("outputs/tables/author_tableS6_APLP1_fraction_author_subtypes_wide.csv"))

  conditions_3 <- data.table::data.table(
    condition_plot = c("non-AD", "AD-AT8+", "AD-AT8-"),
    pct_col = c("pct.2 (non-AD)", "pct.1 (AD-AT8+)", "pct.3 (AD-AT8-)"),
    count_col = c("n_nonAD", "n_AD_AT8pos", "n_AD_AT8neg")
  )

  fig3_dt <- data.table::rbindlist(lapply(row_order, function(display_subtype) {
    source_subtypes <- if (display_subtype == "Ex4&5") c("Ex4", "Ex5") else display_subtype
    data.table::rbindlist(lapply(seq_len(nrow(conditions_3)), function(i) {
      res <- get_weighted_fraction(source_subtypes, conditions_3$pct_col[i], conditions_3$count_col[i])
      data.table::data.table(
        neuronal_subtype = display_subtype,
        source_subtypes = paste(source_subtypes, collapse = "+"),
        condition = conditions_3$condition_plot[i],
        APLP1_positive_fraction = res$fraction,
        author_n_cells = res$n_cells,
        source_column = conditions_3$pct_col[i],
        combination_method = res$method
      )
    }))
  }))

  fig3_dt[, neuronal_subtype := factor(neuronal_subtype, levels = rev(row_order))]
  fig3_dt[, condition := factor(condition, levels = c("non-AD", "AD-AT8+", "AD-AT8-"))]
  write_table(fig3_dt, project_path("outputs/tables/author_tableS6_APLP1_fraction_author_subtypes_3groups.csv"))

  wide3_dt <- data.table::dcast(
    fig3_dt,
    neuronal_subtype + source_subtypes ~ condition,
    value.var = "APLP1_positive_fraction"
  )
  wide3_dt[, neuronal_subtype := as.character(neuronal_subtype)]
  wide3_dt[, row_order_i := match(neuronal_subtype, row_order)]
  data.table::setorder(wide3_dt, row_order_i)
  wide3_dt[, row_order_i := NULL]
  write_table(wide3_dt, project_path("outputs/tables/author_tableS6_APLP1_fraction_author_subtypes_3groups_wide.csv"))

  fdr_star <- function(p) {
    if (is.na(p)) {
      ""
    } else if (p < 0.001) {
      "***"
    } else if (p < 0.01) {
      "**"
    } else if (p < 0.05) {
      "*"
    } else {
      ""
    }
  }
  primary_fdr <- data.table::rbindlist(lapply(row_order, function(display_subtype) {
    source_subtypes <- if (display_subtype == "Ex4&5") c("Ex4", "Ex5") else display_subtype
    p_values <- if (primary_p_col %in% names(aplp1_rows)) {
      numeric_clean(aplp1_rows[subtype %in% source_subtypes, get(primary_p_col)])
    } else {
      NA_real_
    }
    p_values <- p_values[!is.na(p_values)]
    p_adj <- if (length(p_values)) min(p_values) else NA_real_
    data.table::data.table(
      neuronal_subtype = display_subtype,
      source_subtypes = paste(source_subtypes, collapse = "+"),
      AD_AT8pos_vs_AD_AT8neg_FDR = p_adj,
      FDR_star = fdr_star(p_adj),
      FDR_source = if (is.na(p_adj)) {
        "not available in author Table S6 APLP1 row"
      } else if (length(source_subtypes) > 1) {
        "minimum component adjusted p value from author Table S6"
      } else {
        "author Table S6 adjusted p value"
      }
    )
  }))
  primary_fdr[, neuronal_subtype := factor(neuronal_subtype, levels = rev(row_order))]
  write_table(primary_fdr, project_path("outputs/tables/author_tableS6_APLP1_primary_FDR_stars.csv"))

  title_detectable <- "Fraction of neuronal somas with detectable APLP1 transcripts\nin the Otero-Garcia BA9 single-soma dataset"

  p <- ggplot2::ggplot(
    fig_dt,
    ggplot2::aes(x = condition, y = neuronal_subtype, fill = APLP1_positive_fraction)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35, width = 0.96, height = 0.96) +
    ggplot2::scale_fill_gradientn(
      colors = c("#f7f4ff", "#c5baff", "#7567ff", "#2b00d4"),
      limits = c(0.28, 0.70),
      breaks = c(0.30, 0.40, 0.50, 0.60, 0.70),
      labels = scales::number_format(accuracy = 0.1),
      oob = scales::squish,
      name = NULL
    ) +
    ggplot2::scale_x_discrete(labels = c("AD-AT8+", "non-AD\ncontrol")) +
    ggplot2::coord_fixed(ratio = 0.30) +
    ggplot2::labs(
      title = title_detectable,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 7, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 7.4, lineheight = 0.95, color = "black"),
      axis.text.x = ggplot2::element_text(size = 7, color = "black", margin = ggplot2::margin(t = 2)),
      axis.text.y = ggplot2::element_text(size = 7, color = "black", face = "bold"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.25),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.25),
      legend.key.height = grid::unit(0.36, "cm"),
      legend.key.width = grid::unit(0.22, "cm"),
      legend.text = ggplot2::element_text(size = 6.5, color = "black"),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  out_base <- "Fig_APLP1_author_subtype_fraction_heatmap"
  ggplot2::ggsave(
    project_path("outputs/figures", paste0(out_base, ".pdf")),
    p,
    width = 3.25,
    height = 2.45,
    useDingbats = FALSE
  )
  ggplot2::ggsave(
    project_path("outputs/figures", paste0(out_base, ".png")),
    p,
    width = 3.25,
    height = 2.45,
    dpi = 600
  )

  na_tiles <- fig3_dt[is.na(APLP1_positive_fraction)]
  star_tiles <- primary_fdr[FDR_star != ""]
  fig3_dt[, n_cell_label_color := ifelse(!is.na(author_n_cells) & author_n_cells >= 2500, "white", "black")]

  p3 <- ggplot2::ggplot(
    fig3_dt,
    ggplot2::aes(x = condition, y = neuronal_subtype, fill = APLP1_positive_fraction)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35, width = 0.96, height = 0.96) +
    ggplot2::geom_text(
      data = na_tiles,
      ggplot2::aes(label = "n/a"),
      color = "#555555",
      size = 2.0,
      inherit.aes = TRUE
    ) +
    ggplot2::geom_text(
      data = star_tiles,
      ggplot2::aes(x = "AD-AT8+", y = neuronal_subtype, label = FDR_star),
      color = "white",
      size = 2.8,
      fontface = "bold",
      vjust = 0.72,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c("#f7f4ff", "#c5baff", "#7567ff", "#2b00d4"),
      limits = c(0.28, 0.70),
      breaks = c(0.30, 0.40, 0.50, 0.60, 0.70),
      labels = scales::number_format(accuracy = 0.1),
      oob = scales::squish,
      na.value = "#f1f1f1",
      name = NULL
    ) +
    ggplot2::scale_x_discrete(labels = c("non-AD", "AD-AT8+", "AD-AT8-")) +
    ggplot2::coord_fixed(ratio = 0.30) +
    ggplot2::labs(
      title = title_detectable,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 7, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 7.4, lineheight = 0.95, color = "black"),
      axis.text.x = ggplot2::element_text(size = 7, color = "black", margin = ggplot2::margin(t = 2)),
      axis.text.y = ggplot2::element_text(size = 7, color = "black", face = "bold"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.25),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.25),
      legend.key.height = grid::unit(0.36, "cm"),
      legend.key.width = grid::unit(0.22, "cm"),
      legend.text = ggplot2::element_text(size = 6.5, color = "black"),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  out_base_3 <- "Fig_APLP1_author_subtype_fraction_heatmap_3groups"
  ggplot2::ggsave(
    project_path("outputs/figures", paste0(out_base_3, ".pdf")),
    p3,
    width = 3.65,
    height = 2.45,
    useDingbats = FALSE
  )
  ggplot2::ggsave(
    project_path("outputs/figures", paste0(out_base_3, ".png")),
    p3,
    width = 3.65,
    height = 2.45,
    dpi = 600
  )

  p3_fraction_panel <- p3 +
    ggplot2::labs(title = "Detectable APLP1 fraction")

  p3_ncells <- ggplot2::ggplot(
    fig3_dt,
    ggplot2::aes(x = condition, y = neuronal_subtype, fill = author_n_cells)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35, width = 0.96, height = 0.96) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(is.na(author_n_cells), "", scales::comma(author_n_cells, accuracy = 1)),
        color = n_cell_label_color
      ),
      size = 1.85,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
      trans = "log10",
      labels = scales::comma,
      name = "n_cells"
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_x_discrete(labels = c("non-AD", "AD-AT8+", "AD-AT8-")) +
    ggplot2::coord_fixed(ratio = 0.30) +
    ggplot2::labs(title = "Author subtype n_cells", x = NULL, y = NULL) +
    ggplot2::theme_classic(base_size = 7, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 7.4, color = "black"),
      axis.text.x = ggplot2::element_text(size = 7, color = "black", margin = ggplot2::margin(t = 2)),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_line(color = "black", linewidth = 0.25),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.25),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )

  p3_combo <- (p3_fraction_panel | p3_ncells) +
    patchwork::plot_layout(widths = c(1.2, 1.0)) +
    patchwork::plot_annotation(
      title = title_detectable,
      caption = "Fraction = raw APLP1 UMI counts > 0. Stars mark author Table S6 adjusted p values for AD-AT8+ vs AD-AT8-: * FDR < 0.05, ** < 0.01, *** < 0.001. Grey/n/a = not available in author Table S6."
    ) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 8.2, color = "black", lineheight = 0.95),
      plot.caption = ggplot2::element_text(size = 5.7, color = "#333333", hjust = 0)
    )

  out_base_3_combo <- "Fig_APLP1_author_subtype_fraction_ncells_3groups"
  ggplot2::ggsave(
    project_path("outputs/figures", paste0(out_base_3_combo, ".pdf")),
    p3_combo,
    width = 7.0,
    height = 2.85,
    useDingbats = FALSE
  )
  ggplot2::ggsave(
    project_path("outputs/figures", paste0(out_base_3_combo, ".png")),
    p3_combo,
    width = 7.0,
    height = 2.85,
    dpi = 600
  )

  legend_text <- c(
    "Figure. APLP1-positive fraction using Otero-Garcia et al. author neuronal subtype labels.",
    "",
    "The heatmaps show the author Table S6 pct values for cells with detectable APLP1 transcripts, defined as raw APLP1 UMI counts > 0, in non-AD, AD-AT8+, and AD-AT8- neuronal somas. Ex4&5 is the cell-count weighted mean of Ex4 and Ex5 using the author Table S6 summary soma counts. The public author supplementary workbooks Table S2 and Table S6 were inspected for barcode-level subtype metadata; no barcode-to-subtype column or sheet was detected. Therefore these figures reproduce the author subtype-level APLP1 fractions from Table S6 rather than reassigning individual GEO H5 barcodes. Stars mark author Table S6 adjusted p values for AD-AT8+ vs AD-AT8- (* FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001). Missing grey tiles indicate values not available in the author Table S6 APLP1 row."
  )
  writeLines(legend_text, project_path("outputs/report/Fig_APLP1_author_subtype_legend.txt"))

  search_note <- c(
    "Author barcode-level subtype metadata search",
    "",
    paste0("Inspected: ", basename(table_s2), " and ", basename(table_s6), "."),
    "Result: no barcode-level or barcode-to-subtype metadata table was detected in these public supplementary workbooks.",
    "Usable author labels: Table S6 provides subtype-level sheets and APLP1 pct columns for AD-AT8+, non-AD, and AD-AT8- groups.",
    "Figure source: author Table S6 APLP1 pct columns; Ex4&5 is cell-count weighted from Ex4 and Ex5 using Table S6 summary counts.",
    "Labeling rule: the control column is written as non-AD control, not as an AT8-negative group."
  )
  writeLines(search_note, file.path(metadata_dir, "metadata_search_summary.txt"))

  log_message("Author subtype APLP1 figure written to outputs/figures.", log = "10_make_author_subtype_aplp1_figure.log")
})

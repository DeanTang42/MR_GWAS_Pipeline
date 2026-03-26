#!/usr/bin/env Rscript

suppressWarnings({
  parse_args <- function(args) {
    parsed <- list()
    i <- 1
    while (i <= length(args)) {
      key <- args[[i]]
      if (!startsWith(key, "--")) {
        stop("Unexpected argument: ", key, call. = FALSE)
      }
      if (i == length(args)) {
        stop("Missing value for argument: ", key, call. = FALSE)
      }
      parsed[[substring(key, 3)]] <- args[[i + 1]]
      i <- i + 2
    }
    parsed
  }

  args <- parse_args(commandArgs(trailingOnly = TRUE))
  required <- c(
    "r-lib-path", "input", "output", "mr-role",
    "clump-r2", "clump-kb", "clump-p1", "clump-pop",
    "plink", "bfile"
  )
  missing <- required[!required %in% names(args)]
  if (length(missing) > 0) {
    stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  .libPaths(c(args[["r-lib-path"]], .libPaths()))

  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(TwoSampleMR))

  raw <- fread(args[["input"]])
  raw$bim_id <- toupper(raw$bim_id)
  raw$variant_id <- toupper(raw$variant_id)

  if ("phenotype" %in% names(args)) {
    raw$phenotype <- args[["phenotype"]]
  }
  if ("sample-size" %in% names(args)) {
    raw$N <- as.numeric(args[["sample-size"]])
  }

  format_args <- list(
    dat = as.data.frame(raw),
    type = args[["mr-role"]],
    snp_col = "bim_id",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
  )

  if ("phenotype" %in% names(raw)) {
    format_args$phenotype_col <- "phenotype"
  }
  if ("N" %in% names(raw)) {
    format_args$samplesize_col <- "N"
  }

  formatted <- do.call(format_data, format_args)
  formatted$SNP <- toupper(formatted$SNP)

  if (args[["mr-role"]] == "exposure") {
    clumped <- tryCatch(
      clump_data(
        formatted,
        clump_r2 = as.numeric(args[["clump-r2"]]),
        clump_kb = as.numeric(args[["clump-kb"]]),
        clump_p1 = as.numeric(args[["clump-p1"]]),
        pop = args[["clump-pop"]],
        plink = args[["plink"]],
        bfile = args[["bfile"]]
      ),
      error = function(e) {
        warning("clump_data failed: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(clumped) || nrow(clumped) == 0) {
      formatted <- formatted[0, ]
    } else {
      clumped$SNP <- toupper(clumped$SNP)
      formatted <- formatted[formatted$SNP %in% unique(clumped$SNP), , drop = FALSE]
    }
  }

  annotation_cols <- c("bim_id", "variant_id")
  if ("rsid" %in% names(raw)) {
    annotation_cols <- c(annotation_cols, "rsid")
  }
  annotation <- unique(raw[, ..annotation_cols])
  final <- merge(
    as.data.table(formatted),
    annotation,
    by.x = "SNP",
    by.y = "bim_id",
    all.x = TRUE,
    sort = FALSE
  )

  fwrite(final, args[["output"]])

  cat("MR role:", args[["mr-role"]], "\n")
  cat("Output :", args[["output"]], "\n")
  cat("Rows   :", nrow(final), "\n")
})

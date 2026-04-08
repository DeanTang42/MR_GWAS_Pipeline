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

  get_arg <- function(name, default = NULL) {
    if (name %in% names(args) && nzchar(args[[name]])) {
      return(args[[name]])
    }
    default
  }

  normalize_role <- function(role) {
    role <- tolower(role)
    if (role %in% c("exp", "exposure")) {
      return("exposure")
    }
    if (role %in% c("out", "outcome")) {
      return("outcome")
    }
    stop("Unsupported mr-role: ", role, call. = FALSE)
  }

  fill_eaf <- function(eaf, maf, use_maf_as_eaf) {
    if (!use_maf_as_eaf) {
      return(eaf)
    }
    if (is.null(maf)) {
      return(eaf)
    }
    if (is.null(eaf)) {
      return(maf)
    }
    eaf[is.na(eaf)] <- maf[is.na(eaf)]
    eaf
  }

  .libPaths(c(args[["r-lib-path"]], .libPaths()))

  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(TwoSampleMR))

  role <- normalize_role(args[["mr-role"]])
  input_format <- tolower(get_arg("input-format", "raw"))
  use_maf_as_eaf <- tolower(get_arg("use-maf-as-eaf", "false")) %in% c("1", "true", "yes", "y")
  source_dat <- fread(args[["input"]])

  if (input_format == "standardized") {
    required_cols <- c(
      "BIM_ID", "VARIANT_ID", "CHR", "POS",
      "EFFECT_ALLELE", "OTHER_ALLELE", "BETA", "SE", "P"
    )
    missing_cols <- required_cols[!required_cols %in% names(source_dat)]
    if (length(missing_cols) > 0) {
      stop(
        "Standardized input is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    eaf <- if ("EAF" %in% names(source_dat)) as.numeric(source_dat[["EAF"]]) else rep(NA_real_, nrow(source_dat))
    maf <- if ("MAF" %in% names(source_dat)) as.numeric(source_dat[["MAF"]]) else NULL
    eaf <- fill_eaf(eaf, maf, use_maf_as_eaf)

    raw <- data.table(
      bim_id = source_dat[["BIM_ID"]],
      variant_id = source_dat[["VARIANT_ID"]],
      chromosome = source_dat[["CHR"]],
      base_pair_location = source_dat[["POS"]],
      effect_allele = source_dat[["EFFECT_ALLELE"]],
      other_allele = source_dat[["OTHER_ALLELE"]],
      effect_allele_frequency = eaf,
      beta = as.numeric(source_dat[["BETA"]]),
      standard_error = as.numeric(source_dat[["SE"]]),
      p_value = as.numeric(source_dat[["P"]])
    )
    if ("RSID" %in% names(source_dat)) {
      raw$rsid <- source_dat[["RSID"]]
    }
  } else if (input_format == "raw") {
    raw <- source_dat
  } else {
    stop("Unsupported input-format: ", input_format, call. = FALSE)
  }

  if ("region-chr" %in% names(args) || "region-start" %in% names(args) || "region-end" %in% names(args)) {
    region_required <- c("region-chr", "region-start", "region-end")
    region_missing <- region_required[!region_required %in% names(args)]
    if (length(region_missing) > 0) {
      stop("Missing region arguments: ", paste(region_missing, collapse = ", "), call. = FALSE)
    }
    region_chr <- toupper(gsub("^CHR", "", as.character(args[["region-chr"]]), ignore.case = TRUE))
    region_start <- as.numeric(args[["region-start"]])
    region_end <- as.numeric(args[["region-end"]])
    raw <- raw[
      toupper(gsub("^CHR", "", as.character(chromosome), ignore.case = TRUE)) == region_chr &
        as.numeric(base_pair_location) >= region_start &
        as.numeric(base_pair_location) <= region_end
    ]
  }

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
    type = role,
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

  if (role == "exposure") {
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

  cat("MR role:", role, "\n")
  cat("Input format:", input_format, "\n")
  if ("region-chr" %in% names(args)) {
    cat("Region :", paste0(args[["region-chr"]], ":", args[["region-start"]], "-", args[["region-end"]]), "\n")
  }
  cat("Output :", args[["output"]], "\n")
  cat("Rows   :", nrow(final), "\n")
})

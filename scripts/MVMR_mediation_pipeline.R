#!/usr/bin/env Rscript
# ==============================================================================
# X-M-Y 中介 MVMR 分析流水线
# 用法: Rscript scripts/MVMR_mediation_pipeline.R --x-exp-name X --m-exp-name M \
#         --x-std X.tsv.gz --m-std M.tsv.gz --y-std Y.tsv.gz [其他参数]
# ==============================================================================

get_script_path <- function() {
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg) == 0) {
    return(normalizePath(getwd(), mustWork = FALSE))
  }
  normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)
}

parse_cli_args <- function(args) {
  parsed <- list(
    project_dir = normalizePath(file.path(dirname(get_script_path()), ".."), mustWork = FALSE),
    r_lib_path = Sys.getenv("MR_PIPELINE_R_LIB_PATH", unset = "/home/ding/R/4.4.1_MR"),
    x_exp_name = NULL,
    m_exp_name = NULL,
    x_std = NULL,
    m_std = NULL,
    y_std = NULL,
    x_label = NULL,
    m_label = NULL,
    y_label = NULL,
    output_dir = NULL,
    help = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("-h", "--help")) {
      parsed$help <- TRUE
      i <- i + 1
    } else if (arg == "--project-dir") {
      parsed$project_dir <- args[[i + 1]]
      i <- i + 2
    } else if (arg == "--r-lib-path") {
      parsed$r_lib_path <- args[[i + 1]]
      i <- i + 2
    } else if (arg == "--x-exp-name") {
      parsed$x_exp_name <- args[[i + 1]]
      i <- i + 2
    } else if (arg == "--m-exp-name") {
      parsed$m_exp_name <- args[[i + 1]]
      i <- i + 2
    } else if (arg == "--x-std") {
      parsed$x_std <- args[[i + 1]]
      i <- i + 2
    } else if (arg == "--m-std") {
      parsed$m_std <- args[[i + 1]]
      i <- i + 2
    } else if (arg == "--y-std") {
      parsed$y_std <- args[[i + 1]]
      i <- i + 2
    } else if (arg == "--x-label") {
      parsed$x_label <- args[[i + 1]]
      i <- i + 2
    } else if (arg == "--m-label") {
      parsed$m_label <- args[[i + 1]]
      i <- i + 2
    } else if (arg == "--y-label") {
      parsed$y_label <- args[[i + 1]]
      i <- i + 2
    } else if (arg == "--output-dir") {
      parsed$output_dir <- args[[i + 1]]
      i <- i + 2
    } else if (startsWith(arg, "--")) {
      stop("未知参数: ", arg, call. = FALSE)
    } else {
      stop("未识别的位置参数: ", arg, call. = FALSE)
    }
  }

  parsed
}

print_usage <- function() {
  cat(
    "用法:\n",
    "  Rscript scripts/MVMR_mediation_pipeline.R \\\n",
    "    --x-exp-name <X_exp_name> --m-exp-name <M_exp_name> \\\n",
    "    --x-std <X_standardized.tsv.gz> --m-std <M_standardized.tsv.gz> --y-std <Y_standardized.tsv.gz> \\\n",
    "    --x-label <X名称> --m-label <M名称> --y-label <Y名称> \\\n",
    "    --output-dir <结果目录>\n\n",
    "说明:\n",
    "  1. X/M 的 IV 清单来自 data/exp 下的 exposure 文件。\n",
    "  2. X/M/Y 的 beta/se/alleles/eaf 从 standardized 文件中提取。\n",
    "  3. X-Y 与 M-Y 使用 TwoSampleMR::harmonise_data() 配对协调后，构建 MVMR 宽表。\n",
    sep = ""
  )
}

cli <- parse_cli_args(commandArgs(trailingOnly = TRUE))
required_args <- c(
  cli$x_exp_name, cli$m_exp_name, cli$x_std, cli$m_std, cli$y_std,
  cli$x_label, cli$m_label, cli$y_label, cli$output_dir
)
if (cli$help || any(vapply(required_args, is.null, logical(1)))) {
  print_usage()
  quit(status = if (cli$help) 0 else 1)
}

.libPaths(c(cli$r_lib_path, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(MVMR)
})

project_dir <- normalizePath(cli$project_dir, mustWork = FALSE)
exposure_dir <- Sys.getenv(
  "MR_PIPELINE_EXP_DIR",
  unset = Sys.getenv("MR_PIPELINE_EXPOSURE_DIR", unset = file.path(project_dir, "data", "exp"))
)
result_dir <- normalizePath(cli$output_dir, mustWork = FALSE)
clump_plink <- Sys.getenv("MR_PIPELINE_MVMR_CLUMP_PLINK", unset = "/home/ding/miniconda3/envs/GWAS/bin/plink")
clump_bfile <- Sys.getenv("MR_PIPELINE_MVMR_CLUMP_BFILE", unset = "/home/ding/MR_LPA/Ref/g1000_eur/g1000_eur_colon")
clump_r2 <- as.numeric(Sys.getenv("MR_PIPELINE_MVMR_CLUMP_R2", unset = "0.01"))
clump_kb <- as.numeric(Sys.getenv("MR_PIPELINE_MVMR_CLUMP_KB", unset = "1000"))
clump_p1 <- as.numeric(Sys.getenv("MR_PIPELINE_MVMR_CLUMP_P1", unset = "1"))

safe_num_scalar <- function(x, index = 1) {
  if (is.null(x) || length(x) < index) {
    return(NA_real_)
  }
  value <- suppressWarnings(as.numeric(x[[index]]))
  if (length(value) == 0) {
    return(NA_real_)
  }
  value[[1]]
}

first_non_missing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

write_txt <- function(df, path) {
  fwrite(df, file = path, sep = "\t", na = "NA", quote = FALSE)
}

parse_clumped_snps <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    return(character())
  }
  clumped <- fread(path, fill = TRUE)
  if (!"SNP" %in% names(clumped)) {
    return(character())
  }
  unique(toupper(as.character(clumped$SNP)))
}

run_global_clump <- function(input_dt, out_dir, prefix = "global_clump") {
  clump_input_path <- file.path(out_dir, paste0(prefix, "_input.txt"))
  clump_out_prefix <- file.path(out_dir, prefix)
  write_txt(input_dt, clump_input_path)

  cmd <- c(
    "--bfile", clump_bfile,
    "--clump", clump_input_path,
    "--clump-snp-field", "SNP",
    "--clump-field", "P",
    "--clump-kb", as.character(clump_kb),
    "--clump-r2", as.character(clump_r2),
    "--clump-p1", as.character(clump_p1),
    "--out", clump_out_prefix
  )
  output <- system2(clump_plink, args = cmd, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  writeLines(output, paste0(clump_out_prefix, ".log"))
  if (status != 0L) {
    stop("PLINK global clump 失败，请检查: ", paste0(clump_out_prefix, ".log"), call. = FALSE)
  }
  lead_snps <- parse_clumped_snps(paste0(clump_out_prefix, ".clumped"))
  if (length(lead_snps) == 0) {
    stop("全局 clump 未返回 lead SNP，请检查: ", paste0(clump_out_prefix, ".clumped"), call. = FALSE)
  }
  lead_snps
}

fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

fmt_sci <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "e"))
}

resolve_input_file <- function(dir_path, dataset_name) {
  candidates <- c(
    file.path(dir_path, paste0(dataset_name, ".csv")),
    file.path(dir_path, paste0(dataset_name, ".csv.gz"))
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("未找到 exposure 文件: ", paste(candidates, collapse = " / "), call. = FALSE)
  }
  normalizePath(existing[[1]], mustWork = TRUE)
}

normalize_exp_dat <- function(dt, label_fallback) {
  required <- c(
    "SNP", "beta.exposure", "se.exposure",
    "effect_allele.exposure", "other_allele.exposure", "pval.exposure"
  )
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("exposure 文件缺少必要列: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  dt <- copy(dt)
  dt[, SNP := toupper(as.character(SNP))]
  dt[, effect_allele.exposure := toupper(as.character(effect_allele.exposure))]
  dt[, other_allele.exposure := toupper(as.character(other_allele.exposure))]
  if ("eaf.exposure" %in% names(dt)) {
    dt[, eaf.exposure := suppressWarnings(as.numeric(eaf.exposure))]
  } else {
    dt[, eaf.exposure := NA_real_]
  }
  if (!"exposure" %in% names(dt)) {
    dt[, exposure := label_fallback]
  }
  if (!"id.exposure" %in% names(dt)) {
    dt[, id.exposure := label_fallback]
  }
  dt
}

load_standardized_subset <- function(path, snps) {
  available <- names(fread(path, nrows = 0))
  required <- c("BIM_ID", "VARIANT_ID", "RSID", "EFFECT_ALLELE", "OTHER_ALLELE", "BETA", "SE", "P", "EAF")
  select_cols <- intersect(required, available)
  missing <- setdiff(c("BIM_ID", "EFFECT_ALLELE", "OTHER_ALLELE", "BETA", "SE", "P"), select_cols)
  if (length(missing) > 0) {
    stop("standardized 文件缺少必要列: ", paste(missing, collapse = ", "), " -> ", path, call. = FALSE)
  }
  dt <- fread(path, select = select_cols, na.strings = c("NA", ""))
  setDT(dt)
  dt[, BIM_ID := toupper(as.character(BIM_ID))]
  dt <- dt[BIM_ID %in% snps]
  if (!"VARIANT_ID" %in% names(dt)) dt[, VARIANT_ID := NA_character_]
  if (!"RSID" %in% names(dt)) dt[, RSID := NA_character_]
  if (!"EAF" %in% names(dt)) dt[, EAF := NA_real_]
  dt[, EFFECT_ALLELE := toupper(as.character(EFFECT_ALLELE))]
  dt[, OTHER_ALLELE := toupper(as.character(OTHER_ALLELE))]
  dt[, BETA := suppressWarnings(as.numeric(BETA))]
  dt[, SE := suppressWarnings(as.numeric(SE))]
  dt[, P := suppressWarnings(as.numeric(P))]
  dt[, EAF := suppressWarnings(as.numeric(EAF))]
  unique(dt, by = "BIM_ID")
}

to_exposure_std <- function(dt, label) {
  data.table(
    SNP = dt$BIM_ID,
    beta.exposure = dt$BETA,
    se.exposure = dt$SE,
    effect_allele.exposure = dt$EFFECT_ALLELE,
    other_allele.exposure = dt$OTHER_ALLELE,
    eaf.exposure = dt$EAF,
    pval.exposure = dt$P,
    exposure = label,
    id.exposure = label
  )
}

to_outcome_std <- function(dt, label) {
  data.table(
    SNP = dt$BIM_ID,
    beta.outcome = dt$BETA,
    se.outcome = dt$SE,
    effect_allele.outcome = dt$EFFECT_ALLELE,
    other_allele.outcome = dt$OTHER_ALLELE,
    eaf.outcome = dt$EAF,
    pval.outcome = dt$P,
    outcome = label,
    id.outcome = label
  )
}

harmonise_pair <- function(exp_dat, out_dat, pair_label) {
  harm <- harmonise_data(exp_dat, out_dat, action = 2)
  harm <- as.data.table(harm)
  keep <- if ("mr_keep" %in% names(harm)) harm$mr_keep else rep(TRUE, nrow(harm))
  keep[is.na(keep)] <- FALSE
  kept <- harm[keep]
  list(full = harm, keep = kept, label = pair_label)
}

run_ivw_pair <- function(exp_dat, out_dat, pair_label) {
  harm_res <- harmonise_pair(exp_dat, out_dat, pair_label)
  if (nrow(harm_res$keep) == 0) {
    return(list(
      harmonised = harm_res,
      beta = NA_real_,
      se = NA_real_,
      pval = NA_real_,
      nsnp = 0
    ))
  }
  mr_res <- mr(harm_res$keep, method_list = c("mr_ivw"))
  list(
    harmonised = harm_res,
    beta = safe_num_scalar(mr_res$b),
    se = safe_num_scalar(mr_res$se),
    pval = safe_num_scalar(mr_res$pval),
    nsnp = nrow(harm_res$keep)
  )
}

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
cat("输出目录:", result_dir, "\n")

x_exp_file <- resolve_input_file(exposure_dir, cli$x_exp_name)
m_exp_file <- resolve_input_file(exposure_dir, cli$m_exp_name)

cat("正在读取 exposure IV 文件...\n")
x_exp <- normalize_exp_dat(fread(x_exp_file), cli$x_label)
m_exp <- normalize_exp_dat(fread(m_exp_file), cli$m_label)

x_iv <- unique(x_exp$SNP)
m_iv <- unique(m_exp$SNP)
iv_union <- data.table(SNP = sort(unique(c(x_iv, m_iv))))
iv_union[, in_x := SNP %in% x_iv]
iv_union[, in_m := SNP %in% m_iv]
iv_union[, iv_source := fifelse(in_x & in_m, "both", fifelse(in_x, "X", "M"))]
cat("  X IV 数:", length(x_iv), "\n")
cat("  M IV 数:", length(m_iv), "\n")
cat("  IV 并集数:", nrow(iv_union), "\n")

cat("正在从 standardized 文件提取 X/M/Y 数据...\n")
x_std <- load_standardized_subset(cli$x_std, iv_union$SNP)
m_std <- load_standardized_subset(cli$m_std, iv_union$SNP)
y_std <- load_standardized_subset(cli$y_std, iv_union$SNP)
cat("  X standardized 命中:", nrow(x_std), "\n")
cat("  M standardized 命中:", nrow(m_std), "\n")
cat("  Y standardized 命中:", nrow(y_std), "\n")

write_txt(iv_union, file.path(result_dir, "iv_union.txt"))

meta_all <- rbindlist(list(
  x_std[, .(SNP = BIM_ID, variant_id = VARIANT_ID, rsid = RSID)],
  m_std[, .(SNP = BIM_ID, variant_id = VARIANT_ID, rsid = RSID)],
  y_std[, .(SNP = BIM_ID, variant_id = VARIANT_ID, rsid = RSID)]
), use.names = TRUE, fill = TRUE)
meta_dt <- meta_all[, .(
  variant_id = first_non_missing(variant_id),
  rsid = first_non_missing(rsid)
), by = SNP]

cat("正在进行 X-Y / M-Y 配对协调，用于 MVMR...\n")
xy_pair <- harmonise_pair(to_exposure_std(x_std, cli$x_label), to_outcome_std(y_std, cli$y_label), "X-Y")
my_pair <- harmonise_pair(to_exposure_std(m_std, cli$m_label), to_outcome_std(y_std, cli$y_label), "M-Y")

write_txt(xy_pair$full, file.path(result_dir, "harmonised_xy_full.txt"))
write_txt(xy_pair$keep, file.path(result_dir, "harmonised_xy_keep.txt"))
write_txt(my_pair$full, file.path(result_dir, "harmonised_my_full.txt"))
write_txt(my_pair$keep, file.path(result_dir, "harmonised_my_keep.txt"))

xy_mvmr <- xy_pair$keep[, .(
  SNP,
  betaX1 = beta.exposure,
  sebetaX1 = se.exposure,
  pX1 = pval.exposure,
  betaYG = beta.outcome,
  sebetaYG = se.outcome,
  effect_allele_y_xy = effect_allele.outcome,
  other_allele_y_xy = other_allele.outcome
)]
my_mvmr <- my_pair$keep[, .(
  SNP,
  betaX2 = beta.exposure,
  sebetaX2 = se.exposure,
  pX2 = pval.exposure,
  effect_allele_y_my = effect_allele.outcome,
  other_allele_y_my = other_allele.outcome
)]

mvmr_input <- merge(xy_mvmr, my_mvmr, by = "SNP")

# harmonise_data() 会把 outcome 对齐到各自 exposure 的 effect allele。
# 因此 X-Y 与 M-Y 两次结果里，Y 的等位基因方向可能“相同”或“相反”。
# 对于相反的情况，不应直接丢弃，而应把 M 侧 beta 翻转到 X-Y 的方向。
mvmr_input[, orientation_relation := fifelse(
  effect_allele_y_xy == effect_allele_y_my &
    other_allele_y_xy == other_allele_y_my,
  "same",
  fifelse(
    effect_allele_y_xy == other_allele_y_my &
      other_allele_y_xy == effect_allele_y_my,
    "opposite",
    "mismatch"
  )
)]

mvmr_input <- mvmr_input[orientation_relation != "mismatch"]

if (nrow(mvmr_input) == 0) {
  stop("X-Y 与 M-Y 协调后没有共同且方向一致的 SNP，无法进行 MVMR。", call. = FALSE)
}

orientation_same_n <- sum(mvmr_input$orientation_relation == "same", na.rm = TRUE)
orientation_opposite_n <- sum(mvmr_input$orientation_relation == "opposite", na.rm = TRUE)

mvmr_input[orientation_relation == "opposite", betaX2 := -betaX2]

mvmr_input <- merge(mvmr_input, iv_union[, .(SNP, iv_source)], by = "SNP", all.x = TRUE)
mvmr_input <- merge(mvmr_input, meta_dt, by = "SNP", all.x = TRUE)

p_matrix <- as.matrix(mvmr_input[, .(pX1, pX2)])
storage.mode(p_matrix) <- "numeric"
p_matrix[is.na(p_matrix)] <- Inf
p_clump <- apply(p_matrix, 1, min)
p_clump[!is.finite(p_clump)] <- 1
mvmr_input[, p_clump := p_clump]

setcolorder(mvmr_input, c(
  "SNP", "variant_id", "rsid", "iv_source", "p_clump",
  "betaYG", "sebetaYG", "betaX1", "betaX2", "sebetaX1", "sebetaX2",
  "pX1", "pX2", "effect_allele_y_xy", "other_allele_y_xy", "orientation_relation"
))

write_txt(mvmr_input, file.path(result_dir, "pre_clump_mvmr_input.txt"))
pre_clump_n <- nrow(mvmr_input)
lead_snps <- run_global_clump(mvmr_input[, .(SNP, P = p_clump)], result_dir, prefix = "global_clump")
write_txt(merge(data.table(SNP = sort(unique(lead_snps))), mvmr_input[, .(SNP, variant_id, rsid, iv_source, p_clump)], by = "SNP", all.x = TRUE), file.path(result_dir, "global_clump_leads.txt"))
write_txt(mvmr_input[!SNP %in% lead_snps, .(SNP, variant_id, rsid, iv_source, p_clump)], file.path(result_dir, "global_clump_removed.txt"))
mvmr_input <- mvmr_input[SNP %in% lead_snps]
post_clump_n <- nrow(mvmr_input)
write_txt(mvmr_input, file.path(result_dir, "post_clump_mvmr_input.txt"))
write_txt(mvmr_input, file.path(result_dir, "mvmr_input.txt"))
cat("  协调后待全局 clump SNP 数:", pre_clump_n, "\n")
cat("  全局 clump 后保留 SNP 数:", post_clump_n, "\n")

cat("正在进行 MVMR 分析...\n")
formatted_mvmr <- format_mvmr(
  BXGs = as.data.frame(mvmr_input[, .(betaX1, betaX2)]),
  BYG = mvmr_input$betaYG,
  seBXGs = as.data.frame(mvmr_input[, .(sebetaX1, sebetaX2)]),
  seBYG = mvmr_input$sebetaYG,
  RSID = mvmr_input$SNP
)

strength_res <- strength_mvmr(r_input = formatted_mvmr, gencov = 0)
pleiotropy_res <- pleiotropy_mvmr(r_input = formatted_mvmr, gencov = 0)
ivw_mvmr_res <- ivw_mvmr(r_input = formatted_mvmr, gencov = 0)

writeLines(capture.output(print(strength_res)), file.path(result_dir, "mvmr_strength.txt"))
writeLines(capture.output(print(pleiotropy_res)), file.path(result_dir, "mvmr_pleiotropy.txt"))
writeLines(capture.output(print(ivw_mvmr_res)), file.path(result_dir, "mvmr_ivw_raw.txt"))

ivw_df <- as.data.table(as.data.frame(ivw_mvmr_res), keep.rownames = "exposure_row")
if (!"Estimate" %in% names(ivw_df)) {
  stop("ivw_mvmr() 返回结构异常，未找到 Estimate 列。", call. = FALSE)
}
ivw_df[, exposure_label := fifelse(exposure_row == "exposure1", cli$x_label,
                            fifelse(exposure_row == "exposure2", cli$m_label, exposure_row))]
write_txt(ivw_df, file.path(result_dir, "mvmr_ivw.txt"))

direct_x_beta <- safe_num_scalar(ivw_df[exposure_row == "exposure1", Estimate])
direct_x_se <- safe_num_scalar(ivw_df[exposure_row == "exposure1", `Std. Error`])
direct_x_p <- safe_num_scalar(ivw_df[exposure_row == "exposure1", `Pr(>|t|)`])
direct_m_beta <- safe_num_scalar(ivw_df[exposure_row == "exposure2", Estimate])
direct_m_se <- safe_num_scalar(ivw_df[exposure_row == "exposure2", `Std. Error`])
direct_m_p <- safe_num_scalar(ivw_df[exposure_row == "exposure2", `Pr(>|t|)`])

cat("正在计算 X->M 与 X->Y 总效应...\n")
y_std_for_x <- y_std[BIM_ID %in% x_exp$SNP]
m_std_for_x <- m_std[BIM_ID %in% x_exp$SNP]

x_to_y_total <- run_ivw_pair(x_exp, to_outcome_std(y_std_for_x, cli$y_label), "X->Y")
x_to_m_total <- run_ivw_pair(x_exp, to_outcome_std(m_std_for_x, cli$m_label), "X->M")

write_txt(x_to_y_total$harmonised$full, file.path(result_dir, "harmonised_total_xy_full.txt"))
write_txt(x_to_y_total$harmonised$keep, file.path(result_dir, "harmonised_total_xy_keep.txt"))
write_txt(x_to_m_total$harmonised$full, file.path(result_dir, "harmonised_total_xm_full.txt"))
write_txt(x_to_m_total$harmonised$keep, file.path(result_dir, "harmonised_total_xm_keep.txt"))

indirect_beta <- x_to_m_total$beta * direct_m_beta
indirect_se <- if (is.na(x_to_m_total$se) || is.na(direct_m_se)) {
  NA_real_
} else {
  sqrt((direct_m_beta^2) * (x_to_m_total$se^2) + (x_to_m_total$beta^2) * (direct_m_se^2))
}
indirect_z <- if (is.na(indirect_se) || indirect_se == 0) NA_real_ else indirect_beta / indirect_se
indirect_p <- if (is.na(indirect_z)) NA_real_ else 2 * pnorm(-abs(indirect_z))
prop_mediated <- if (is.na(x_to_y_total$beta) || x_to_y_total$beta == 0) NA_real_ else indirect_beta / x_to_y_total$beta

mediation_summary <- data.table(
  effect = c(
    paste0(cli$x_label, " -> ", cli$y_label, " total"),
    paste0(cli$x_label, " -> ", cli$m_label, " total"),
    paste0(cli$x_label, " -> ", cli$y_label, " direct | ", cli$m_label),
    paste0(cli$m_label, " -> ", cli$y_label, " direct | ", cli$x_label),
    paste0(cli$x_label, " -> ", cli$y_label, " indirect via ", cli$m_label)
  ),
  estimate = c(
    x_to_y_total$beta,
    x_to_m_total$beta,
    direct_x_beta,
    direct_m_beta,
    indirect_beta
  ),
  se = c(
    x_to_y_total$se,
    x_to_m_total$se,
    direct_x_se,
    direct_m_se,
    indirect_se
  ),
  pval = c(
    x_to_y_total$pval,
    x_to_m_total$pval,
    direct_x_p,
    direct_m_p,
    indirect_p
  ),
  nsnp = c(
    x_to_y_total$nsnp,
    x_to_m_total$nsnp,
    nrow(mvmr_input),
    nrow(mvmr_input),
    nrow(mvmr_input)
  )
)
write_txt(mediation_summary, file.path(result_dir, "mediation_summary.txt"))

strength_lines <- capture.output(print(strength_res))
pleiotropy_lines <- capture.output(print(pleiotropy_res))

report_path <- file.path(result_dir, "report.txt")
report_lines <- c(
  "================================================================================",
  paste0("X-M-Y 中介 MVMR 报告: ", cli$x_label, " -> ", cli$m_label, " -> ", cli$y_label),
  "================================================================================",
  "",
  "一、输入与 IV 构建",
  "--------------------------------------------------------------------------------",
  paste0("  X exposure 文件: ", x_exp_file),
  paste0("  M exposure 文件: ", m_exp_file),
  paste0("  X standardized 文件: ", normalizePath(cli$x_std, mustWork = TRUE)),
  paste0("  M standardized 文件: ", normalizePath(cli$m_std, mustWork = TRUE)),
  paste0("  Y standardized 文件: ", normalizePath(cli$y_std, mustWork = TRUE)),
  paste0("  X IV 数: ", length(x_iv)),
  paste0("  M IV 数: ", length(m_iv)),
  paste0("  IV 并集数: ", nrow(iv_union)),
  "",
  "二、配对协调与 MVMR 输入",
  "--------------------------------------------------------------------------------",
  paste0("  X-Y 协调后可用 SNP 数: ", nrow(xy_pair$keep)),
  paste0("  M-Y 协调后可用 SNP 数: ", nrow(my_pair$keep)),
  paste0("  两次协调方向相同 SNP 数: ", orientation_same_n),
  paste0("  两次协调方向相反并已翻转保留 SNP 数: ", orientation_opposite_n),
  paste0("  协调后待全局 clump SNP 数: ", pre_clump_n),
  paste0("  全局 clump 参数: r2=", clump_r2, ", kb=", clump_kb, ", p1=", clump_p1),
  paste0("  全局 clump 后最终 SNP 数: ", post_clump_n),
  "",
  "三、条件工具强度",
  "--------------------------------------------------------------------------------",
  strength_lines,
  "",
  "四、横向多效性检验",
  "--------------------------------------------------------------------------------",
  pleiotropy_lines,
  "",
  "五、MVMR 直接效应",
  "--------------------------------------------------------------------------------",
  paste0("  ", cli$x_label, " -> ", cli$y_label, " | ", cli$m_label, ": beta = ", fmt_num(direct_x_beta, 4),
         ", se = ", fmt_num(direct_x_se, 4), ", P = ", fmt_sci(direct_x_p, 2)),
  paste0("  ", cli$m_label, " -> ", cli$y_label, " | ", cli$x_label, ": beta = ", fmt_num(direct_m_beta, 4),
         ", se = ", fmt_num(direct_m_se, 4), ", P = ", fmt_sci(direct_m_p, 2)),
  "",
  "六、中介分解",
  "--------------------------------------------------------------------------------",
  paste0("  ", cli$x_label, " -> ", cli$y_label, " 总效应: beta = ", fmt_num(x_to_y_total$beta, 4),
         ", se = ", fmt_num(x_to_y_total$se, 4), ", P = ", fmt_sci(x_to_y_total$pval, 2)),
  paste0("  ", cli$x_label, " -> ", cli$m_label, " 总效应: beta = ", fmt_num(x_to_m_total$beta, 4),
         ", se = ", fmt_num(x_to_m_total$se, 4), ", P = ", fmt_sci(x_to_m_total$pval, 2)),
  paste0("  经 ", cli$m_label, " 的间接效应: beta = ", fmt_num(indirect_beta, 4),
         ", se = ", fmt_num(indirect_se, 4), ", P = ", fmt_sci(indirect_p, 2)),
  paste0("  中介比例(间接/总效应): ", fmt_num(prop_mediated, 4)),
  "",
  "七、输出文件",
  "--------------------------------------------------------------------------------",
  paste0("  IV 并集: ", file.path(result_dir, "iv_union.txt")),
  paste0("  全局 clump 输入: ", file.path(result_dir, "global_clump_input.txt")),
  paste0("  全局 clump lead SNP: ", file.path(result_dir, "global_clump_leads.txt")),
  paste0("  全局 clump 剔除 SNP: ", file.path(result_dir, "global_clump_removed.txt")),
  paste0("  clump 前 MVMR 输入: ", file.path(result_dir, "pre_clump_mvmr_input.txt")),
  paste0("  MVMR 输入: ", file.path(result_dir, "mvmr_input.txt")),
  paste0("  MVMR IVW: ", file.path(result_dir, "mvmr_ivw.txt")),
  paste0("  MVMR 工具强度: ", file.path(result_dir, "mvmr_strength.txt")),
  paste0("  MVMR 多效性: ", file.path(result_dir, "mvmr_pleiotropy.txt")),
  paste0("  中介汇总: ", file.path(result_dir, "mediation_summary.txt")),
  "================================================================================"
)
writeLines(report_lines, report_path)

cat("分析完成!\n")
cat("报告文件:", report_path, "\n")

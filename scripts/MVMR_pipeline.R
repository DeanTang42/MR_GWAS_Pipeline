#!/usr/bin/env Rscript
# ==============================================================================
# 通用多 exposure MVMR 分析流水线
# 用法:
#   Rscript scripts/MVMR_pipeline.R \
#     --exp-name exp1 --exp-name exp2 \
#     --exp-std exp1.tsv.gz --exp-std exp2.tsv.gz \
#     --exp-label EXP1 --exp-label EXP2 \
#     --y-std outcome.tsv.gz --y-label OUTCOME \
#     --output-dir results/...
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
    exp_names = character(),
    exp_stds = character(),
    exp_labels = character(),
    y_std = NULL,
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
    } else if (arg == "--exp-name") {
      parsed$exp_names <- c(parsed$exp_names, args[[i + 1]])
      i <- i + 2
    } else if (arg == "--exp-std") {
      parsed$exp_stds <- c(parsed$exp_stds, args[[i + 1]])
      i <- i + 2
    } else if (arg == "--exp-label") {
      parsed$exp_labels <- c(parsed$exp_labels, args[[i + 1]])
      i <- i + 2
    } else if (arg == "--y-std") {
      parsed$y_std <- args[[i + 1]]
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
    "  Rscript scripts/MVMR_pipeline.R \\\n",
    "    --exp-name exp1 --exp-name exp2 [--exp-name exp3 ...] \\\n",
    "    --exp-std exp1.tsv.gz --exp-std exp2.tsv.gz [--exp-std exp3.tsv.gz ...] \\\n",
    "    --exp-label EXP1 --exp-label EXP2 [--exp-label EXP3 ...] \\\n",
    "    --y-std outcome.tsv.gz --y-label OUTCOME \\\n",
    "    --output-dir <结果目录>\n\n",
    "说明:\n",
    "  1. exposure 的 IV 清单来自 data/exp 下对应的 exposure 文件。\n",
    "  2. 所有 exposure 与 outcome 的 beta/se/alleles/eaf 从 standardized 文件中提取。\n",
    "  3. 每个 exposure 都与 outcome 做一次 TwoSampleMR::harmonise_data(action = 2)。\n",
    "  4. 以第 1 个 exposure 为锚点，若其他 exposure 的 outcome 方向相反，则翻转相应 beta 后保留。\n",
    sep = ""
  )
}

cli <- parse_cli_args(commandArgs(trailingOnly = TRUE))
if (cli$help) {
  print_usage()
  quit(status = 0)
}

if (length(cli$exp_names) < 2) {
  print_usage()
  stop("至少需要 2 个 exposure。", call. = FALSE)
}
if (length(cli$exp_names) != length(cli$exp_stds)) {
  stop("--exp-name 与 --exp-std 数量不一致。", call. = FALSE)
}
if (length(cli$exp_labels) == 0) {
  cli$exp_labels <- cli$exp_names
}
if (length(cli$exp_names) != length(cli$exp_labels)) {
  stop("--exp-name 与 --exp-label 数量不一致。", call. = FALSE)
}
if (is.null(cli$y_std) || is.null(cli$y_label) || is.null(cli$output_dir)) {
  print_usage()
  stop("缺少必要参数: --y-std / --y-label / --output-dir", call. = FALSE)
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

slugify <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
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

harmonise_pair <- function(exp_dat, out_dat) {
  harm <- harmonise_data(exp_dat, out_dat, action = 2)
  harm <- as.data.table(harm)
  keep <- if ("mr_keep" %in% names(harm)) harm$mr_keep else rep(TRUE, nrow(harm))
  keep[is.na(keep)] <- FALSE
  list(full = harm, keep = harm[keep])
}

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
cat("输出目录:", result_dir, "\n")

exp_n <- length(cli$exp_names)
exp_files <- vapply(cli$exp_names, function(name) resolve_input_file(exposure_dir, name), character(1))

cat("正在读取 exposure IV 文件...\n")
exp_dat_list <- vector("list", exp_n)
for (i in seq_len(exp_n)) {
  exp_dat_list[[i]] <- normalize_exp_dat(fread(exp_files[[i]]), cli$exp_labels[[i]])
}

iv_sources_long <- rbindlist(lapply(seq_len(exp_n), function(i) {
  data.table(SNP = unique(exp_dat_list[[i]]$SNP), exposure_label = cli$exp_labels[[i]])
}))
iv_union <- unique(iv_sources_long[, .(SNP)])
setorder(iv_union, SNP)
iv_sources <- iv_sources_long[, .(iv_sources = paste(sort(unique(exposure_label)), collapse = ",")), by = SNP]

cat("  Exposure 数量:", exp_n, "\n")
for (i in seq_len(exp_n)) {
  cat("   -", cli$exp_labels[[i]], "原始 IV 数:", length(unique(exp_dat_list[[i]]$SNP)), "\n")
}
cat("  IV 并集数:", nrow(iv_union), "\n")
write_txt(merge(iv_union, iv_sources, by = "SNP", all.x = TRUE), file.path(result_dir, "iv_union.txt"))

cat("正在从 standardized 文件提取各 exposure 与 outcome 数据...\n")
exp_std_list <- vector("list", exp_n)
for (i in seq_len(exp_n)) {
  exp_std_list[[i]] <- load_standardized_subset(cli$exp_stds[[i]], iv_union$SNP)
  cat("   -", cli$exp_labels[[i]], "standardized 命中:", nrow(exp_std_list[[i]]), "\n")
}
y_std <- load_standardized_subset(cli$y_std, iv_union$SNP)
cat("   -", cli$y_label, "standardized 命中:", nrow(y_std), "\n")

meta_all <- rbindlist(c(
  lapply(exp_std_list, function(dt) dt[, .(SNP = BIM_ID, variant_id = VARIANT_ID, rsid = RSID)]),
  list(y_std[, .(SNP = BIM_ID, variant_id = VARIANT_ID, rsid = RSID)])
), use.names = TRUE, fill = TRUE)
meta_dt <- meta_all[, .(
  variant_id = first_non_missing(variant_id),
  rsid = first_non_missing(rsid)
), by = SNP]

cat("正在进行 pairwise harmonise...\n")
pair_results <- vector("list", exp_n)
pair_counts <- vector("list", exp_n)
for (i in seq_len(exp_n)) {
  pair_results[[i]] <- harmonise_pair(to_exposure_std(exp_std_list[[i]], cli$exp_labels[[i]]), to_outcome_std(y_std, cli$y_label))
  slug <- sprintf("%02d_%s", i, slugify(cli$exp_labels[[i]]))
  write_txt(pair_results[[i]]$full, file.path(result_dir, paste0("harmonised_", slug, "_full.txt")))
  write_txt(pair_results[[i]]$keep, file.path(result_dir, paste0("harmonised_", slug, "_keep.txt")))
  pair_counts[[i]] <- data.table(
    exposure_index = i,
    exposure_label = cli$exp_labels[[i]],
    original_iv_n = length(unique(exp_dat_list[[i]]$SNP)),
    harmonised_all_n = nrow(pair_results[[i]]$full),
    harmonised_keep_n = nrow(pair_results[[i]]$keep)
  )
}
pair_counts_dt <- rbindlist(pair_counts, use.names = TRUE, fill = TRUE)
write_txt(pair_counts_dt, file.path(result_dir, "pairwise_harmonise_summary.txt"))

common_snps <- Reduce(intersect, lapply(pair_results, function(item) item$keep$SNP))
if (length(common_snps) == 0) {
  stop("所有 exposure 与 outcome 协调后没有共同 SNP，无法进行 MVMR。", call. = FALSE)
}

baseline <- pair_results[[1]]$keep[SNP %in% common_snps, .(
  SNP,
  betaYG = beta.outcome,
  sebetaYG = se.outcome,
  effect_allele_y_ref = effect_allele.outcome,
  other_allele_y_ref = other_allele.outcome
)]

mvmr_input <- baseline
orientation_summary <- list()

for (i in seq_len(exp_n)) {
  cur <- pair_results[[i]]$keep[SNP %in% common_snps, .(
    SNP,
    betaX = beta.exposure,
    sebetaX = se.exposure,
    pX = pval.exposure,
    effect_allele_y_cur = effect_allele.outcome,
    other_allele_y_cur = other_allele.outcome
  )]
  mvmr_input <- merge(mvmr_input, cur, by = "SNP")
  relation_col <- paste0("orientation_relation_", i)
  mvmr_input[, (relation_col) := fifelse(
    effect_allele_y_ref == effect_allele_y_cur &
      other_allele_y_ref == other_allele_y_cur,
    "same",
    fifelse(
      effect_allele_y_ref == other_allele_y_cur &
        other_allele_y_ref == effect_allele_y_cur,
      "opposite",
      "mismatch"
    )
  )]
  if (i > 1) {
    mvmr_input[get(relation_col) == "opposite", betaX := -betaX]
  }
  orientation_summary[[i]] <- data.table(
    exposure_index = i,
    exposure_label = cli$exp_labels[[i]],
    orientation_same_n = sum(mvmr_input[[relation_col]] == "same", na.rm = TRUE),
    orientation_opposite_n = sum(mvmr_input[[relation_col]] == "opposite", na.rm = TRUE),
    orientation_mismatch_n = sum(mvmr_input[[relation_col]] == "mismatch", na.rm = TRUE)
  )
  setnames(mvmr_input, c("betaX", "sebetaX", "pX"), c(paste0("betaX", i), paste0("sebetaX", i), paste0("pX", i)))
  mvmr_input[, c("effect_allele_y_cur", "other_allele_y_cur") := NULL]
}

orientation_dt <- rbindlist(orientation_summary, use.names = TRUE, fill = TRUE)
write_txt(orientation_dt, file.path(result_dir, "orientation_summary.txt"))

mismatch_cols <- grep("^orientation_relation_", names(mvmr_input), value = TRUE)
if (length(mismatch_cols) > 0) {
  mismatch_dt <- mvmr_input[, ..mismatch_cols]
  keep_mask <- rowSums(as.data.frame(mismatch_dt) == "mismatch", na.rm = TRUE) == 0
  mvmr_input <- mvmr_input[keep_mask]
}
if (nrow(mvmr_input) == 0) {
  stop("共同 SNP 在方向统一后全部失效，无法进行 MVMR。", call. = FALSE)
}

mvmr_input <- merge(mvmr_input, iv_sources, by = "SNP", all.x = TRUE)
mvmr_input <- merge(mvmr_input, meta_dt, by = "SNP", all.x = TRUE)

p_cols <- paste0("pX", seq_len(exp_n))
p_matrix <- as.matrix(mvmr_input[, ..p_cols])
storage.mode(p_matrix) <- "numeric"
p_matrix[is.na(p_matrix)] <- Inf
p_clump <- apply(p_matrix, 1, min)
p_clump[!is.finite(p_clump)] <- 1
mvmr_input[, p_clump := p_clump]

setcolorder(mvmr_input, c("SNP", "variant_id", "rsid", "iv_sources", "p_clump", setdiff(names(mvmr_input), c("SNP", "variant_id", "rsid", "iv_sources", "p_clump"))))
write_txt(mvmr_input, file.path(result_dir, "pre_clump_mvmr_input.txt"))

pre_clump_n <- nrow(mvmr_input)
lead_snps <- run_global_clump(mvmr_input[, .(SNP, P = p_clump)], result_dir, prefix = "global_clump")
write_txt(merge(data.table(SNP = sort(unique(lead_snps))), mvmr_input[, .(SNP, variant_id, rsid, iv_sources, p_clump)], by = "SNP", all.x = TRUE), file.path(result_dir, "global_clump_leads.txt"))
write_txt(mvmr_input[!SNP %in% lead_snps, .(SNP, variant_id, rsid, iv_sources, p_clump)], file.path(result_dir, "global_clump_removed.txt"))
mvmr_input <- mvmr_input[SNP %in% lead_snps]
post_clump_n <- nrow(mvmr_input)
write_txt(mvmr_input, file.path(result_dir, "post_clump_mvmr_input.txt"))
write_txt(mvmr_input, file.path(result_dir, "mvmr_input.txt"))
cat("  协调后待全局 clump SNP 数:", pre_clump_n, "\n")
cat("  全局 clump 后保留 SNP 数:", post_clump_n, "\n")

bx_cols <- paste0("betaX", seq_len(exp_n))
sebx_cols <- paste0("sebetaX", seq_len(exp_n))
formatted_mvmr <- format_mvmr(
  BXGs = as.data.frame(mvmr_input[, ..bx_cols]),
  BYG = mvmr_input$betaYG,
  seBXGs = as.data.frame(mvmr_input[, ..sebx_cols]),
  seBYG = mvmr_input$sebetaYG,
  RSID = mvmr_input$SNP
)

cat("正在进行 MVMR 分析...\n")
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
row_map <- data.table(exposure_row = paste0("exposure", seq_len(exp_n)), exposure_label = cli$exp_labels)
ivw_df <- merge(ivw_df, row_map, by = "exposure_row", all.x = TRUE)
write_txt(ivw_df, file.path(result_dir, "mvmr_ivw.txt"))

report_lines <- c(
  "================================================================================",
  paste0("通用 MVMR 报告: ", paste(cli$exp_labels, collapse = " + "), " -> ", cli$y_label),
  "================================================================================",
  "",
  "一、输入与 IV 构建",
  "--------------------------------------------------------------------------------",
  paste0("  Outcome standardized 文件: ", normalizePath(cli$y_std, mustWork = TRUE)),
  paste0("  Outcome 名称: ", cli$y_label),
  paste0("  Exposure 数量: ", exp_n),
  vapply(seq_len(exp_n), function(i) {
    paste0("  Exposure ", i, ": ", cli$exp_labels[[i]],
           " | exp 文件=", exp_files[[i]],
           " | standardized 文件=", normalizePath(cli$exp_stds[[i]], mustWork = TRUE))
  }, character(1)),
  vapply(seq_len(exp_n), function(i) {
    paste0("    原始 IV 数(", cli$exp_labels[[i]], "): ", length(unique(exp_dat_list[[i]]$SNP)))
  }, character(1)),
  paste0("  IV 并集数: ", nrow(iv_union)),
  "",
  "二、Pairwise 协调",
  "--------------------------------------------------------------------------------",
  vapply(seq_len(exp_n), function(i) {
    paste0("  ", cli$exp_labels[[i]], " -> ", cli$y_label,
           " 协调后可用 SNP 数: ", nrow(pair_results[[i]]$keep))
  }, character(1)),
  "",
  "三、方向统一与最终输入",
  "--------------------------------------------------------------------------------",
  vapply(seq_len(exp_n), function(i) {
    paste0("  ", cli$exp_labels[[i]],
           " 相对第1个 exposure 的方向: same=",
           orientation_dt[exposure_index == i, orientation_same_n],
           ", opposite=", orientation_dt[exposure_index == i, orientation_opposite_n],
           ", mismatch=", orientation_dt[exposure_index == i, orientation_mismatch_n])
  }, character(1)),
  paste0("  协调后待全局 clump SNP 数: ", pre_clump_n),
  paste0("  全局 clump 参数: r2=", clump_r2, ", kb=", clump_kb, ", p1=", clump_p1),
  paste0("  全局 clump 后最终 SNP 数: ", post_clump_n),
  "",
  "四、条件工具强度",
  "--------------------------------------------------------------------------------",
  capture.output(print(strength_res)),
  "",
  "五、横向多效性检验",
  "--------------------------------------------------------------------------------",
  capture.output(print(pleiotropy_res)),
  "",
  "六、MVMR 直接效应",
  "--------------------------------------------------------------------------------",
  vapply(seq_len(nrow(ivw_df)), function(i) {
    paste0("  ", ivw_df$exposure_label[[i]],
           " -> ", cli$y_label,
           ": beta = ", fmt_num(safe_num_scalar(ivw_df$Estimate[[i]]), 4),
           ", se = ", fmt_num(safe_num_scalar(ivw_df$`Std. Error`[[i]]), 4),
           ", P = ", fmt_sci(safe_num_scalar(ivw_df$`Pr(>|t|)`[[i]]), 2))
  }, character(1)),
  "",
  "七、输出文件",
  "--------------------------------------------------------------------------------",
  paste0("  IV 并集: ", file.path(result_dir, "iv_union.txt")),
  paste0("  Pairwise 协调汇总: ", file.path(result_dir, "pairwise_harmonise_summary.txt")),
  paste0("  方向汇总: ", file.path(result_dir, "orientation_summary.txt")),
  paste0("  全局 clump 输入: ", file.path(result_dir, "global_clump_input.txt")),
  paste0("  全局 clump lead SNP: ", file.path(result_dir, "global_clump_leads.txt")),
  paste0("  全局 clump 剔除 SNP: ", file.path(result_dir, "global_clump_removed.txt")),
  paste0("  clump 前 MVMR 输入: ", file.path(result_dir, "pre_clump_mvmr_input.txt")),
  paste0("  MVMR 输入: ", file.path(result_dir, "mvmr_input.txt")),
  paste0("  MVMR IVW: ", file.path(result_dir, "mvmr_ivw.txt")),
  paste0("  MVMR 工具强度: ", file.path(result_dir, "mvmr_strength.txt")),
  paste0("  MVMR 多效性: ", file.path(result_dir, "mvmr_pleiotropy.txt")),
  "================================================================================"
)
writeLines(report_lines, file.path(result_dir, "report.txt"))

cat("分析完成!\n")
cat("报告文件:", file.path(result_dir, "report.txt"), "\n")

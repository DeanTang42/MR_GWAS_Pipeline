#!/usr/bin/env Rscript
# ==============================================================================
# MR分析流水线 - 孟德尔随机化单暴露单结局分析
# 用法: Rscript MR_single_pipeline.R <暴露名> <结局名>
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
    use_steiger = FALSE,
    steiger_pval = 0.05,
    positional = character(0),
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
    } else if (arg == "--use-steiger") {
      parsed$use_steiger <- TRUE
      i <- i + 1
    } else if (arg == "--steiger-pval") {
      parsed$steiger_pval <- as.numeric(args[[i + 1]])
      i <- i + 2
    } else if (startsWith(arg, "--")) {
      stop("未知参数: ", arg, call. = FALSE)
    } else {
      parsed$positional <- c(parsed$positional, arg)
      i <- i + 1
    }
  }

  parsed
}

print_usage <- function() {
  cat(
    "用法:\n",
    "  Rscript scripts/MR_single_pipeline.R [选项] <暴露名> <结局名>\n\n",
    "可选参数:\n",
    "  --project-dir <dir>   项目根目录，默认自动推断\n",
    "  --r-lib-path <dir>    R 包库路径，默认读取 MR_PIPELINE_R_LIB_PATH 或 /home/ding/R/4.4.1_MR\n",
    "  --use-steiger         启用 Steiger 方向性过滤\n",
    "  --steiger-pval <num>  Steiger 过滤阈值，默认 0.05\n",
    "  -h, --help            显示帮助\n\n",
    "说明:\n",
    "  暴露/结局名需与 data/exposure 和 data/outcome 目录下的 .csv 或 .csv.gz 文件名一致(不含后缀)\n",
    sep = ""
  )
}

cli <- parse_cli_args(commandArgs(trailingOnly = TRUE))
if (cli$help || length(cli$positional) < 2) {
  print_usage()
  quit(status = if (cli$help) 0 else 1)
}

.libPaths(c(cli$r_lib_path, .libPaths()))

# ============================
# 加载必要的R包
# ============================
suppressPackageStartupMessages({
  library(data.table)   # 高效数据读取
  library(TwoSampleMR)  # 两样本MR核心包
  library(MRPRESSO)     # MR-PRESSO离群值检测
  library(mr.raps)      # RAPS方法
  library(dplyr)        # 数据处理
  library(ggplot2)      # 绑定ggplot2
})

# ============================
# 项目目录配置
# ============================
project_dir  <- normalizePath(cli$project_dir, mustWork = FALSE)
exposure_dir <- Sys.getenv("MR_PIPELINE_EXPOSURE_DIR", unset = file.path(project_dir, "data", "exposure"))
outcome_dir  <- Sys.getenv("MR_PIPELINE_OUTCOME_DIR", unset = file.path(project_dir, "data", "outcome"))
result_dir   <- Sys.getenv("MR_PIPELINE_RESULTS_DIR", unset = file.path(project_dir, "results"))

# ============================
# 分析选项配置
# ============================
use_steiger   <- cli$use_steiger
steiger_pval  <- cli$steiger_pval
exp_name <- cli$positional[1]  # 暴露名称
out_name <- cli$positional[2]  # 结局名称

# ============================
# 辅助函数
# ============================

# 安全获取列(列不存在时返回NA)
safe_col <- function(df, name) {
  if (name %in% names(df)) df[[name]] else rep(NA, nrow(df))
}

# 计算F统计量
calc_fstat <- function(beta, se) {
  ifelse(is.na(beta) | is.na(se) | se == 0, NA, (beta^2) / (se^2))
}

# 数值格式化
fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

# 科学计数法格式化
fmt_sci <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "e"))
}

# 写入制表符分隔文件
write_txt <- function(df, path) {
  fwrite(df, file = path, sep = "\t", na = "NA", quote = FALSE)
}

resolve_input_file <- function(dir_path, dataset_name) {
  candidates <- c(
    file.path(dir_path, paste0(dataset_name, ".csv")),
    file.path(dir_path, paste0(dataset_name, ".csv.gz"))
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop(
      "未找到输入文件: ",
      paste(candidates, collapse = " / "),
      call. = FALSE
    )
  }
  normalizePath(existing[1], mustWork = TRUE)
}

# ============================
# 主分析流程
# ============================
run_mr_analysis <- function(exp_name, out_name) {
  
  # ---------------------------
  # 1. 创建输出目录
  # ---------------------------
  output_dir <- file.path(result_dir, paste0(exp_name, "_", out_name))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat("输出目录:", output_dir, "\n")
  
  # ---------------------------
  # 2. 读取数据(.csv.gz格式)
  # ---------------------------
  cat("正在读取数据...\n")
  exp_file <- resolve_input_file(exposure_dir, exp_name)
  out_file <- resolve_input_file(outcome_dir, out_name)
  
  # fread可以直接读取.gz压缩文件
  Exp <- fread(exp_file)
  Out <- fread(out_file)
  cat("  暴露SNP数:", nrow(Exp), "\n")
  cat("  结局SNP数:", nrow(Out), "\n")
  
  # ---------------------------
  # 3. 数据协调(Harmonization)
  # ---------------------------
  cat("正在进行数据协调...\n")
  HarmonizedData <- harmonise_data(Exp, Out, action = 2)
  HarmonizedData <- as.data.frame(HarmonizedData)
  
  # 筛选mr_keep=TRUE的SNP
  mr_keep <- if ("mr_keep" %in% names(HarmonizedData)) HarmonizedData$mr_keep else rep(TRUE, nrow(HarmonizedData))
  mr_keep[is.na(mr_keep)] <- FALSE
  removed_harmonise <- HarmonizedData$SNP[!mr_keep]
  HarmonizedData <- HarmonizedData[mr_keep, , drop = FALSE]
  n_harmonised <- nrow(HarmonizedData)
  cat("  协调后SNP数:", n_harmonised, "\n")
  
  # ---------------------------
  # 4. Steiger方向性过滤(可选)
  # ---------------------------
  removed_steiger <- character(0)
  steiger_applied <- FALSE
  
  if (use_steiger) {
    cat("正在进行Steiger方向性过滤...\n")
    steiger_res <- tryCatch(steiger_filtering(HarmonizedData), error = function(e) NULL)
    
    if (!is.null(steiger_res)) {
      steiger_applied <- TRUE
      # 筛选steiger_pval < 阈值的SNP(方向正确)
      keep <- steiger_res$steiger_pval < steiger_pval
      keep[is.na(keep)] <- FALSE
      removed_steiger <- steiger_res$SNP[!keep]
      HarmonizedData <- steiger_res[keep, , drop = FALSE]
      cat("  Steiger过滤移除SNP数:", length(removed_steiger), "\n")
      cat("  过滤后SNP数:", nrow(HarmonizedData), "\n")
    } else {
      cat("  Steiger过滤运行失败，跳过\n")
    }
  }
  
  n_after_steiger <- nrow(HarmonizedData)
  
  # ---------------------------
  # 4. MR-PRESSO离群值检测与去除
  # ---------------------------
  cat("正在进行MR-PRESSO离群值检测...\n")
  
  outlier_snps <- character(0)
  mr_presso_global_p_raw <- NA      # 第一轮(原始)全局P值
  mr_presso_global_p_clean <- NA    # 第二轮(去除离群值后)全局P值
  mr_presso_distortion_p <- NA
  mr_presso_ok <- TRUE
  
  # 第一次运行MR-PRESSO检测离群值
  MrPressoRes <- tryCatch({
    mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = HarmonizedData,
      NbDistribution = 3000,
      SignifThreshold = 0.05
    )
  }, error = function(e) NULL)
  
  if (!is.null(MrPressoRes)) {
    mr_presso_global_p_raw <- MrPressoRes$`MR-PRESSO results`$`Global Test`$Pvalue
    
    # 检查是否有离群值
    outliers <- MrPressoRes$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
    if (!is.null(outliers) && is.numeric(outliers) && length(outliers) > 0) {
      outlier_snps <- HarmonizedData$SNP[outliers]
      cat("  检测到", length(outlier_snps), "个离群值SNP，正在去除...\n")
      
      # 去除离群值
      HarmonizedData <- HarmonizedData[-outliers, , drop = FALSE]
      
      # 重新运行MR-PRESSO获取校正后结果
      MrPressoRes <- tryCatch({
        mr_presso(
          BetaOutcome = "beta.outcome",
          BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome",
          SdExposure = "se.exposure",
          OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE,
          data = HarmonizedData,
          NbDistribution = 3000,
          SignifThreshold = 0.05
        )
      }, error = function(e) NULL)
      
      if (!is.null(MrPressoRes)) {
        # 第二轮全局P值(去除离群值后)
        mr_presso_global_p_clean <- MrPressoRes$`MR-PRESSO results`$`Global Test`$Pvalue
        mr_presso_distortion_p <- MrPressoRes$`MR-PRESSO results`$`Distortion Test`$Pvalue
      }
    } else {
      cat("  未检测到离群值\n")
      # 无离群值时，clean P值等于raw P值
      mr_presso_global_p_clean <- mr_presso_global_p_raw
    }
  } else {
    mr_presso_ok <- FALSE
    cat("  MR-PRESSO运行失败\n")
  }
  
  n_final <- nrow(HarmonizedData)
  cat("  最终用于分析的SNP数:", n_final, "\n")
  
  # ---------------------------
  # 5. MR分析(使用去除离群值后的数据)
  # ---------------------------
  cat("正在进行MR分析...\n")
  
  # 主要MR方法(IVW, Egger, Weighted Median)
  MrRes <- mr(HarmonizedData, method_list = c(
    "mr_ivw",               # 逆方差加权法
    "mr_egger_regression",  # MR-Egger回归
    "mr_weighted_median"    # 加权中位数法
  ))
  
  # 异质性检验
  Het <- mr_heterogeneity(HarmonizedData)
  
  # 多效性检验(Egger截距)
  Ple <- mr_pleiotropy_test(HarmonizedData)
  
  # ---------------------------
  # 6. RAPS分析
  # ---------------------------
  cat("正在进行RAPS分析...\n")
  raps_beta <- NA
  raps_se <- NA
  raps_pval <- NA
  
  raps_result <- tryCatch({
    mr.raps.simple.robust(
      HarmonizedData$beta.exposure,
      HarmonizedData$beta.outcome,
      HarmonizedData$se.exposure,
      HarmonizedData$se.outcome
    )
  }, error = function(e) NULL)
  
  if (!is.null(raps_result)) {
    raps_beta <- raps_result$beta.hat
    raps_se <- raps_result$beta.se
    raps_pval <- raps_result$beta.p.value
  }
  
  # ---------------------------
  # 7. 整理所有结果
  # ---------------------------
  Res <- data.frame(
    method = MrRes$method,
    nsnp = MrRes$nsnp,
    beta = as.numeric(MrRes$b),
    se = as.numeric(MrRes$se),
    pval = as.numeric(MrRes$pval),
    stringsAsFactors = FALSE
  )
  
  # 添加MR-PRESSO结果(校正后)
  if (mr_presso_ok && !is.null(MrPressoRes)) {
    # 优先使用校正后结果
    mr_presso_beta <- ifelse(is.na(MrPressoRes$`Main MR results`$`Causal Estimate`[2]),
                              MrPressoRes$`Main MR results`$`Causal Estimate`[1],
                              MrPressoRes$`Main MR results`$`Causal Estimate`[2])
    mr_presso_se <- ifelse(is.na(MrPressoRes$`Main MR results`$Sd[2]),
                            MrPressoRes$`Main MR results`$Sd[1],
                            MrPressoRes$`Main MR results`$Sd[2])
    mr_presso_pval <- ifelse(is.na(MrPressoRes$`Main MR results`$`P-value`[2]),
                              MrPressoRes$`Main MR results`$`P-value`[1],
                              MrPressoRes$`Main MR results`$`P-value`[2])
    
    Res <- rbind(Res, data.frame(
      method = "MR-PRESSO",
      nsnp = n_final,
      beta = mr_presso_beta,
      se = mr_presso_se,
      pval = mr_presso_pval,
      stringsAsFactors = FALSE
    ))
  }
  
  # 添加RAPS结果
  if (!is.na(raps_beta)) {
    Res <- rbind(Res, data.frame(
      method = "RAPS",
      nsnp = n_final,
      beta = raps_beta,
      se = raps_se,
      pval = raps_pval,
      stringsAsFactors = FALSE
    ))
  }
  
  # 计算效应值和置信区间
  Res$OR <- exp(Res$beta)
  Res$"95%CI_lower" <- exp(Res$beta - 1.96 * Res$se)
  Res$"95%CI_upper" <- exp(Res$beta + 1.96 * Res$se)
  
  # ---------------------------
  # 8. 工具变量信息表(最终使用的IV)
  # ---------------------------
  iv_table <- data.frame(
    SNP = HarmonizedData$SNP,
    effect_allele = safe_col(HarmonizedData, "effect_allele.exposure"),
    other_allele = safe_col(HarmonizedData, "other_allele.exposure"),
    EAF = safe_col(HarmonizedData, "eaf.exposure"),
    beta = safe_col(HarmonizedData, "beta.exposure"),
    se = safe_col(HarmonizedData, "se.exposure"),
    pval = safe_col(HarmonizedData, "pval.exposure"),
    F_stat = calc_fstat(safe_col(HarmonizedData, "beta.exposure"), 
                        safe_col(HarmonizedData, "se.exposure")),
    stringsAsFactors = FALSE
  )
  
  # F统计量汇总
  iv_f <- iv_table$F_stat
  f_mean <- mean(iv_f, na.rm = TRUE)
  f_min <- min(iv_f, na.rm = TRUE)
  f_lt10 <- sum(iv_f < 10, na.rm = TRUE)
  
  # ---------------------------
  # 9. 绘图(关键: 防止Rplots.pdf生成)
  # ---------------------------
  cat("正在生成图形...\n")
  
  # 获取单SNP分析和留一法分析结果
  ResSingle <- mr_singlesnp(HarmonizedData)
  ResLoo <- mr_leaveoneout(HarmonizedData)
  
  # 1) 散点回归图
  plot_scatter <- file.path(output_dir, "scatter.png")
  p_scatter <- mr_scatter_plot(MrRes, HarmonizedData)[[1]]
  ggsave(plot_scatter, plot = p_scatter, width = 8, height = 6, dpi = 300, bg = "white")
  
  # 2) 漏斗图
  plot_funnel <- file.path(output_dir, "funnel.png")
  p_funnel <- mr_funnel_plot(ResSingle)[[1]]
  ggsave(plot_funnel, plot = p_funnel, width = 8, height = 6, dpi = 300, bg = "white")
  
  # 3) 留一法图
  plot_loo <- file.path(output_dir, "leaveoneout.png")
  p_loo <- mr_leaveoneout_plot(ResLoo)[[1]]
  ggsave(plot_loo, plot = p_loo, width = 8, height = 6, dpi = 300, bg = "white")
  
  # 关闭所有图形设备，防止生成Rplots.pdf
  graphics.off()
  
  cat("  已生成3张图形\n")
  
  # ---------------------------
  # 10. 保存结果文件
  # ---------------------------
  cat("正在保存结果...\n")
  
  # MR结果表
  results_path <- file.path(output_dir, "mr_results.txt")
  write_txt(Res, results_path)
  
  # 工具变量表
  iv_path <- file.path(output_dir, "iv_table.txt")
  write_txt(iv_table, iv_path)
  
  # 异质性检验结果
  het_path <- file.path(output_dir, "heterogeneity.txt")
  write_txt(Het, het_path)
  
  # 多效性检验结果
  ple_path <- file.path(output_dir, "pleiotropy.txt")
  write_txt(Ple, ple_path)
  
  # 离群值SNP列表
  if (length(outlier_snps) > 0) {
    outlier_path <- file.path(output_dir, "outliers.txt")
    writeLines(outlier_snps, outlier_path)
  }
  
  # Steiger过滤移除的SNP列表
  if (length(removed_steiger) > 0) {
    steiger_removed_path <- file.path(output_dir, "steiger_removed.txt")
    writeLines(removed_steiger, steiger_removed_path)
  }
  
  # ---------------------------
  # 11. 生成中文报告
  # ---------------------------
  het_ivw <- Het[grep("Inverse variance weighted", Het$method), , drop = FALSE]
  
  report_path <- file.path(output_dir, "report.txt")
  report_lines <- c(
    "================================================================================",
    paste0("孟德尔随机化分析报告: ", exp_name, " → ", out_name),
    "================================================================================",
    "",
    "一、分析概况",
    "--------------------------------------------------------------------------------",
    paste0("  暴露: ", exp_name),
    paste0("  结局: ", out_name),
    paste0("  暴露数据文件: ", exp_file),
    paste0("  结局数据文件: ", out_file),
    paste0("  分析时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "二、数据处理",
    "--------------------------------------------------------------------------------",
    paste0("  原始暴露SNP数: ", nrow(Exp)),
    paste0("  原始结局SNP数: ", nrow(Out)),
    paste0("  协调后SNP数: ", n_harmonised),
    paste0("  协调过程移除的SNP数: ", length(removed_harmonise)),
    "",
    if (use_steiger) "三、Steiger方向性过滤" else NULL,
    if (use_steiger) "--------------------------------------------------------------------------------" else NULL,
    if (use_steiger) paste0("  是否启用: 是") else NULL,
    if (use_steiger) paste0("  过滤阈值: P < ", steiger_pval) else NULL,
    if (use_steiger) paste0("  移除的SNP数: ", length(removed_steiger)) else NULL,
    if (use_steiger && length(removed_steiger) > 0) paste0("  移除的SNP: ", paste(removed_steiger, collapse = ", ")) else NULL,
    if (use_steiger) paste0("  过滤后SNP数: ", n_after_steiger) else NULL,
    if (use_steiger) "" else NULL,
    if (use_steiger) "四、MR-PRESSO离群值检测" else "三、MR-PRESSO离群值检测",
    "--------------------------------------------------------------------------------",
    paste0("  原始全局检验P值(去除前): ", fmt_sci(mr_presso_global_p_raw, 2)),
    paste0("  检测到离群值数量: ", length(outlier_snps)),
    if (length(outlier_snps) > 0) paste0("  离群值SNP: ", paste(outlier_snps, collapse = ", ")) else NULL,
    paste0("  去除离群值后SNP数: ", n_final),
    paste0("  校正后全局检验P值(去除后): ", fmt_sci(mr_presso_global_p_clean, 2)),
    "",
    if (use_steiger) "五、工具变量(IV)信息" else "四、工具变量(IV)信息",
    "--------------------------------------------------------------------------------",
    paste0("  最终用于分析的IV数量: ", n_final),
    paste0("  F统计量均值: ", fmt_num(f_mean, 2)),
    paste0("  F统计量最小值: ", fmt_num(f_min, 2)),
    paste0("  F<10的IV数量: ", f_lt10),
    "",
    if (use_steiger) "六、MR分析结果" else "五、MR分析结果",
    "--------------------------------------------------------------------------------"
  )
  
  # 添加结果表格
  for (i in 1:nrow(Res)) {
    report_lines <- c(report_lines,
      paste0("  ", Res$method[i], ":"),
      paste0("    beta = ", fmt_num(Res$beta[i], 4), 
             ", OR = ", fmt_num(Res$OR[i], 4),
             ", 95%CI = [", fmt_num(Res$`95%CI_lower`[i], 4), ", ", 
             fmt_num(Res$`95%CI_upper`[i], 4), "]"),
      paste0("    P值 = ", fmt_sci(Res$pval[i], 2)),
      ""
    )
  }
  
  report_lines <- c(report_lines,
    if (use_steiger) "七、质量控制检验" else "六、质量控制检验",
    "--------------------------------------------------------------------------------",
    paste0("  异质性检验 (IVW): Q = ", 
           ifelse(nrow(het_ivw) > 0, fmt_num(het_ivw$Q[1], 2), "NA"),
           ", P = ", 
           ifelse(nrow(het_ivw) > 0, fmt_sci(het_ivw$Q_pval[1], 2), "NA")),
    paste0("  多效性检验 (Egger截距): intercept = ", fmt_num(Ple$egger_intercept[1], 4),
           ", P = ", fmt_sci(Ple$pval[1], 2)),
    paste0("  MR-PRESSO全局检验(去除离群值后): P = ", fmt_sci(mr_presso_global_p_clean, 2)),
    if (!is.na(mr_presso_distortion_p)) paste0("  MR-PRESSO扭曲检验: P = ", fmt_sci(mr_presso_distortion_p, 2)) else NULL,
    "",
    if (use_steiger) "八、输出文件" else "七、输出文件",
    "--------------------------------------------------------------------------------",
    paste0("  MR结果表: ", results_path),
    paste0("  工具变量表: ", iv_path),
    paste0("  异质性检验: ", het_path),
    paste0("  多效性检验: ", ple_path),
    if (length(outlier_snps) > 0) paste0("  离群值列表: ", file.path(output_dir, "outliers.txt")) else NULL,
    if (length(removed_steiger) > 0) paste0("  Steiger过滤移除SNP: ", file.path(output_dir, "steiger_removed.txt")) else NULL,
    paste0("  散点回归图: ", plot_scatter),
    paste0("  漏斗图: ", plot_funnel),
    paste0("  留一法图: ", plot_loo),
    "",
    "================================================================================"
  )
  
  # 移除NULL元素
  report_lines <- report_lines[!sapply(report_lines, is.null)]
  writeLines(report_lines, report_path)
  
  cat("\n分析完成!\n")
  cat("结果目录:", output_dir, "\n")
  cat("报告文件:", report_path, "\n")
  
  return(output_dir)
}

# ============================
# 执行分析
# ============================
run_mr_analysis(exp_name, out_name)

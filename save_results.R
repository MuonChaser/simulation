# save_results.R
# 汇总模拟结果并保存
#
# save_results(results, label) -> stats_df
#   results : N_simu × 8 矩阵（来自 run_simu()）
#   label   : 实验标识符，如 "g1_K25"
#
# 输出文件（写入 results/ 目录）:
#   results/{label}_results.rds  — 原始结果矩阵
#   results/{label}_stats.csv    — 汇总统计量
#
# 使用 rsimsum 计算 Bias / SD / RMSE / CP / CI_Length 及其 MC 标准误

save_results <- function(results, label) {
  dir.create("results", showWarnings = FALSE)
  rds_path <- file.path("results", paste0(label, "_results.rds"))
  csv_path <- file.path("results", paste0(label, "_stats.csv"))

  # 保存原始矩阵
  saveRDS(results, rds_path)
  cat("Saved:", rds_path, "\n")

  # 转为长格式（rsimsum 所需：bias 已是 estimate - truth）
  long_df <- rbind(
    data.frame(bias = results[, "bias_strata_a0"], se = results[, "se_strata_a0"],
               Method = "Strata", A = "A=0"),
    data.frame(bias = results[, "bias_group_a0"],  se = results[, "se_group_a0"],
               Method = "Group",  A = "A=0"),
    data.frame(bias = results[, "bias_strata_a1"], se = results[, "se_strata_a1"],
               Method = "Strata", A = "A=1"),
    data.frame(bias = results[, "bias_group_a1"],  se = results[, "se_group_a1"],
               Method = "Group",  A = "A=1")
  )

  # rsimsum: true = 0（bias 已去中心化）
  ss <- simsum(
    data       = long_df,
    estvarname = "bias",
    true       = 0,
    se         = "se",
    methodvar  = "Method",
    ref        = "Group",
    by         = "A"
  )
  ss_tidy <- tidy(ss)

  # 提取统计量并构建 stats_df
  get_stat <- function(stat_name, meth, a_val) {
    ss_tidy[ss_tidy$stat == stat_name & ss_tidy$Method == meth & ss_tidy$A == a_val, ]
  }

  stats_list <- list()
  for (meth in c("Strata", "Group")) {
    for (a_val in c("A=0", "A=1")) {
      b  <- get_stat("bias",    meth, a_val)
      es <- get_stat("empse",   meth, a_val)   # 经验 SD
      ms <- get_stat("mse",     meth, a_val)   # MSE; RMSE via delta method
      cv <- get_stat("cover",   meth, a_val)   # 95% CI 覆盖率
      ml <- get_stat("modelse", meth, a_val)   # 模型 SE 均值

      rmse_val  <- sqrt(ms$est)
      rmse_mcse <- if (rmse_val > 0) ms$mcse / (2 * rmse_val) else 0

      stats_list[[paste(meth, a_val)]] <- data.frame(
        Method    = meth,
        A         = a_val,
        Mean      = b$est,
        Mean_mcse = b$mcse,
        SD        = es$est,
        SD_mcse   = es$mcse,
        RMSE      = rmse_val,
        RMSE_mcse = rmse_mcse,
        CP        = cv$est * 100,          # 比例 -> 百分比
        CP_mcse   = cv$mcse * 100,
        CI_Length = ml$est * 2 * 1.96,    # 均值 CI 长度 = 2 × 1.96 × mean(SE)
        CI_mcse   = ml$mcse * 2 * 1.96
      )
    }
  }

  stats_df <- do.call(rbind, stats_list)
  rownames(stats_df) <- NULL

  write.csv(stats_df, csv_path, row.names = FALSE)
  cat("Saved:", csv_path, "\n")
  print(stats_df)
  cat("Experiment", label, "completed!\n")

  invisible(stats_df)
}

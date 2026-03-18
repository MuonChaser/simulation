# simu.R
# 并行 Monte Carlo 模拟
#
# run_simu(data, cade_a0, cade_a1, N_simu) -> results (N_simu × 8 矩阵)
#
# 每行 8 列:
#   bias_strata_a0, bias_group_a0, bias_strata_a1, bias_group_a1
#   se_strata_a0,   se_group_a0,   se_strata_a1,   se_group_a1
#
# 依赖函数（在 main.R 中已加载）:
#   assign_strata()  — 层次感知分配 (proposed)
#   assign_no_strata() — 聚类整体分配 (benchmark)
#   CADEreg_new()    — 提议的 2SLS 估计量
#   CADEreg()        — 基准估计量（来自 experiment 包）

run_simu <- function(data, cade_a0, cade_a1, N_simu) {
  n_cores   <- min(parallel::detectCores() - 1, 120)
  cl        <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)

  # 进度追踪：每完成一次 replication，worker 向文件追加一个字节
  # 在另一个终端用以下命令监控:
  #   watch -n 2 'printf "%d / N_simu\n" $(wc -c < /tmp/simu_progress.tick)'
  tick_file <- "/tmp/simu_progress.tick"
  file.create(tick_file, showWarnings = FALSE)
  writeLines(character(0), tick_file)   # 清空

  cat(sprintf("Starting %d replications on %d cores ...\n", N_simu, n_cores))
  cat(sprintf("  Progress: watch -n 2 'echo \"$(wc -c < %s) / %d done\"'\n",
              tick_file, N_simu))

  results <- foreach(
    i             = 1:N_simu,
    .combine      = rbind,
    .multicombine = TRUE,
    .maxcombine   = N_simu,
    .packages     = "experiment",
    .export       = c("assign_strata", "assign_no_strata", "CADEreg_new",
                      "tick_file")
  ) %dorng% {
    tryCatch(cat(".", file = tick_file, append = TRUE),
             error = function(e) NULL)
    d_strata <- assign_strata(data)
    d_group  <- assign_no_strata(data)
    r_strata <- CADEreg_new(d_strata)
    r_group  <- CADEreg(d_group)
    c(
      r_strata$CADE0 - cade_a0,   # bias_strata_a0
      r_group$CADE0  - cade_a0,   # bias_group_a0
      r_strata$CADE1 - cade_a1,   # bias_strata_a1
      r_group$CADE1  - cade_a1,   # bias_group_a1
      sqrt(r_strata$var0.reg),    # se_strata_a0
      sqrt(r_group$var0.reg),     # se_group_a0
      sqrt(r_strata$var1.reg),    # se_strata_a1
      sqrt(r_group$var1.reg)      # se_group_a1
    )
  }

  cat("Replications done.\n")

  results <- matrix(results, ncol = 8)
  colnames(results) <- c(
    "bias_strata_a0", "bias_group_a0",
    "bias_strata_a1", "bias_group_a1",
    "se_strata_a0",   "se_group_a0",
    "se_strata_a1",   "se_group_a1"
  )
  results
}

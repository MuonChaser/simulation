# main.R
# CADE 模拟研究入口脚本
#
# 三组实验设计:
#   第一组: K=25/50/100,  n_jk=10,  J=10   (小层数量变化)
#   第二组: K=2,          n_jk=100/200/500, J=10   (层内样本量变化)
#   第三组: K=4,          n_jk=20,  J=20/50/100    (聚类数量变化)
#
# 每组实验运行 N_simu 次 Monte Carlo 重复（论文使用 10,000）
# 结果保存至 results/，论文图表写入 paper/

# ---- 加载函数和包 ----
source("CADEreg_new.R")   # CADEreg_new()
source("assign.R")        # assign_strata(), assign_no_strata()
source("gen_data.R")      # gen_data()
source("simu.R")          # run_simu()
source("save_results.R")  # save_results()

library(experiment)      # CADEreg() —— 基准估计量
library(foreach)
library(doParallel)
library(doRNG)
library(rsimsum)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

N_simu <- 10000   # 论文使用 10,000；调试时可改小
set.seed(2026)

# ---- 单组实验函数 ----
run_experiment <- function(K_val, n_jk_val, J_val, label_val) {
  cat("\n========================================\n")
  cat("Running experiment:", label_val, "\n")
  cat("K =", K_val, ", n_jk =", n_jk_val, ", J =", J_val, "\n")
  cat("========================================\n")

  pop     <- gen_data(J_val, K_val, n_jk_val)
  results <- run_simu(pop$data, pop$cade_a0, pop$cade_a1, N_simu)
  stats   <- save_results(results, label_val)
  cat("Finished:", label_val, "at", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

  stats
}

# ---- 运行三组实验 ----
all_results <- list()

cat("\n\n######## 第一组: 变化 K ########\n")
for (K_val in c(25, 50, 100)) {
  label <- paste0("g1_K", K_val)
  all_results[[label]] <- run_experiment(K_val, 10, 10, label)
}

cat("\n\n######## 第二组: 变化 n_jk ########\n")
for (n_jk_val in c(100, 200, 500)) {
  label <- paste0("g2_n", n_jk_val)
  all_results[[label]] <- run_experiment(2, n_jk_val, 10, label)
}

cat("\n\n######## 第三组: 变化 J ########\n")
for (J_val in c(20, 50, 100)) {
  label <- paste0("g3_J", J_val)
  all_results[[label]] <- run_experiment(4, 20, J_val, label)
}

# ---- 保存汇总结果 ----
saveRDS(all_results, "results/all_experiments_summary.rds")
cat("\n\nAll experiments completed! Summary saved to results/all_experiments_summary.rds\n")

# ---- 生成论文图表 ----
cat("\n\n######## 生成论文图表 ########\n")
source("plot_results.R")    # -> paper/simulation_bias.pdf, paper/simulation_sd.pdf
source("table_results.R")   # -> paper/simulation_table.tex
cat("Figures and table saved to paper/\n")

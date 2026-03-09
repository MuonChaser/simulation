# main.R - 三组实验脚本
# 第一组: K=25,50,100; n_jk=10; J=10
# 第二组: K=2; n_jk=100,200,500; J=10
# 第三组: K=4; n_jk=20; J=25,50,100

source("CADEreg_new.R")
source("assign.R")
library(experiment)
library(ggplot2)
library(foreach)
library(doParallel)
library(doRNG)

N_simu <- 10
# set.seed(123)
set.seed(2026)

# 运行单组实验的函数
run_experiment <- function(K_val, n_jk_val, J_val, label_val) {
  cat("\n========================================\n")
  cat("Running experiment:", label_val, "\n")
  cat("K =", K_val, ", n_jk =", n_jk_val, ", J =", J_val, "\n")
  cat("========================================\n")

  # 设置全局变量
  K <<- K_val
  n_jk <<- n_jk_val
  J <<- J_val
  label <<- label_val

  # 生成数据 -> data, N, n.avg, cade_a0, cade_a1
  source("gen_data.R", local = FALSE)

  # 运行模拟 -> results
  source("simu.R", local = FALSE)

  # 保存结果 -> RDS, CSV, PNG
  source("save_results.R", local = FALSE)

  return(stats_df)
}

# ========================================
# 运行三组实验
# ========================================

all_results <- list()

# 第一组: 变化 K (K=25,50,100; n_jk=10; J=10)
cat("\n\n######## 第一组实验: 变化 K ########\n")
for (K_val in c(25, 50, 100)) {
  all_results[[paste0("g1_K", K_val)]] <- run_experiment(K_val, 10, 10, paste0("g1_K", K_val))
}

# 第二组: 变化 n_jk (K=2; n_jk=100,200,500; J=10)
cat("\n\n######## 第二组实验: 变化 n_jk ########\n")
for (n_jk_val in c(100, 200, 500)) {
  all_results[[paste0("g2_n", n_jk_val)]] <- run_experiment(2, n_jk_val, 10, paste0("g2_n", n_jk_val))
}

# 第三组: 变化 J (K=4; n_jk=20; J=20,50,100)
cat("\n\n######## 第三组实验: 变化 J ########\n")
for (J_val in c(20, 50, 100)) {
  all_results[[paste0("g3_J", J_val)]] <- run_experiment(4, 20, J_val, paste0("g3_J", J_val))
}

# 保存汇总结果
saveRDS(all_results, "all_experiments_summary.rds")
cat("\n\nAll experiments completed! Summary saved to all_experiments_summary.rds\n")

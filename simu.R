# simu.R
# 依赖全局变量: data, J, K, N, cade_a0, cade_a1, N_simu
# 依赖函数: assign_strata, assign_no_strata, CADEreg_new, CADEreg
# 输出全局变量: results (N_simu x 4 矩阵)

n_cores <- min(parallel::detectCores() - 1, 100)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

results <- foreach(
  i = 1:N_simu,
  .combine = rbind,
  .packages = "experiment",
  .export = c("data", "J", "K", "N", "cade_a0", "cade_a1",
              "assign_strata", "assign_no_strata", "CADEreg_new")
) %dorng% {
  d_strata <- assign_strata(data)
  d_group <- assign_no_strata(data)
  r_strata <- CADEreg_new(d_strata)
  r_group <- CADEreg(d_group)
  c(r_strata$CADE0 - cade_a0,
    r_group$CADE0 - cade_a0,
    r_strata$CADE1 - cade_a1,
    r_group$CADE1 - cade_a1)
}

stopCluster(cl)

colnames(results) <- c("bias_strata_a0", "bias_group_a0",
                       "bias_strata_a1", "bias_group_a1")

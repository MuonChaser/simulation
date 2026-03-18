# assign.R
# 随机分配函数
#
# assign_strata(data)    — 层内按比例分配 Z（提议方法）
# assign_no_strata(data) — 聚类整体分配 Z（基准方法）
#
# 依赖全局变量: J, K, N（由 gen_data() 返回的总体参数提供，
#   但函数需在对应环境中调用，见 simu.R 中的 foreach 闭包）
#
# 返回: 添加了 A, Z, D, Y 列的 data.frame

assign_strata <- function(data) {
  J <- nlevels(data$id)
  K <- max(data$B)
  N <- nrow(data)

  # 随机将 J/2 个聚类分配到 A=1，其余 A=0
  A_j    <- sample(c(rep(0, J / 2), rep(1, J / 2)))
  data$A <- A_j[as.integer(data$id)]

  z_prop_high <- 0.7
  z_prop_low  <- 0.3

  # 预计算各聚类的行索引，避免循环内重复 which()
  idx_by_j <- split(seq_len(N), data$id)

  data$Z <- rep(0L, N)
  for (j in 1:J) {
    idx_j <- idx_by_j[[j]]
    B_j   <- data$B[idx_j]
    prop  <- if (A_j[j] == 1) z_prop_high else z_prop_low
    for (k in 1:K) {
      idx_jk  <- idx_j[B_j == k]
      n_treat <- round(prop * length(idx_jk))
      if (n_treat > 0L) data$Z[sample(idx_jk, n_treat)] <- 1L
    }
  }

  # 矩阵列索引替代 4 层嵌套 ifelse
  col_idx <- 1L + data$Z + 2L * data$A
  d_mat   <- cbind(data$d_z0_a0, data$d_z1_a0, data$d_z0_a1, data$d_z1_a1)
  y_mat   <- cbind(data$y_z0_a0, data$y_z1_a0, data$y_z0_a1, data$y_z1_a1)
  rows    <- seq_len(N)
  data$D  <- d_mat[cbind(rows, col_idx)]
  data$Y  <- y_mat[cbind(rows, col_idx)]
  data
}


assign_no_strata <- function(data) {
  J <- nlevels(data$id)
  N <- nrow(data)

  # 随机将 J/2 个聚类分配到 A=1，其余 A=0
  A_j    <- sample(c(rep(0, J / 2), rep(1, J / 2)))
  data$A <- A_j[as.integer(data$id)]

  z_prop_high <- 0.7
  z_prop_low  <- 0.3

  # 预计算各聚类的行索引，避免循环内重复 which()
  idx_by_j <- split(seq_len(N), data$id)

  data$Z <- rep(0L, N)
  for (j in 1:J) {
    idx_j   <- idx_by_j[[j]]
    prop    <- if (A_j[j] == 1) z_prop_high else z_prop_low
    n_treat <- round(prop * length(idx_j))
    data$Z[sample(idx_j, n_treat)] <- 1L
  }

  # 矩阵列索引替代 4 层嵌套 ifelse
  col_idx <- 1L + data$Z + 2L * data$A
  d_mat   <- cbind(data$d_z0_a0, data$d_z1_a0, data$d_z0_a1, data$d_z1_a1)
  y_mat   <- cbind(data$y_z0_a0, data$y_z1_a0, data$y_z0_a1, data$y_z1_a1)
  rows    <- seq_len(N)
  data$D  <- d_mat[cbind(rows, col_idx)]
  data$Y  <- y_mat[cbind(rows, col_idx)]
  data
}

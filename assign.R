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

  data$Z <- rep(0L, N)
  for (j in 1:J) {
    idx_j <- which(data$id == j)
    B_j   <- data$B[idx_j]
    prop  <- ifelse(A_j[j] == 1, z_prop_high, z_prop_low)
    for (k in 1:K) {
      idx_jk  <- idx_j[B_j == k]
      n_treat <- round(prop * length(idx_jk))
      data$Z[sample(idx_jk, n_treat)] <- 1L
    }
  }

  data$D <- ifelse(data$Z == 1 & data$A == 1, data$d_z1_a1,
             ifelse(data$Z == 0 & data$A == 1, data$d_z0_a1,
              ifelse(data$Z == 1 & data$A == 0, data$d_z1_a0,
                                                data$d_z0_a0)))
  data$Y <- ifelse(data$Z == 1 & data$A == 1, data$y_z1_a1,
             ifelse(data$Z == 0 & data$A == 1, data$y_z0_a1,
              ifelse(data$Z == 1 & data$A == 0, data$y_z1_a0,
                                                data$y_z0_a0)))
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

  data$Z <- rep(0L, N)
  for (j in 1:J) {
    idx_j   <- which(data$id == j)
    prop    <- ifelse(A_j[j] == 1, z_prop_high, z_prop_low)
    n_treat <- round(prop * length(idx_j))
    data$Z[sample(idx_j, n_treat)] <- 1L
  }

  data$D <- ifelse(data$Z == 1 & data$A == 1, data$d_z1_a1,
             ifelse(data$Z == 0 & data$A == 1, data$d_z0_a1,
              ifelse(data$Z == 1 & data$A == 0, data$d_z1_a0,
                                                data$d_z0_a0)))
  data$Y <- ifelse(data$Z == 1 & data$A == 1, data$y_z1_a1,
             ifelse(data$Z == 0 & data$A == 1, data$y_z0_a1,
              ifelse(data$Z == 1 & data$A == 0, data$y_z1_a0,
                                                data$y_z0_a0)))
  data
}

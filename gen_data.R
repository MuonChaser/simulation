# gen_data.R
# 生成模拟总体数据
#
# gen_data(J, K, n_jk) -> list(data, N, n.avg, cade_a0, cade_a1)
#   J     : 聚类数
#   K     : 每个聚类内的层数
#   n_jk  : 每层内的个体数

gen_data <- function(J, K, n_jk) {
  n.avg <- n_jk * K
  N     <- J * n.avg

  dat       <- list()
  dat$id    <- as.factor(rep(1:J, each = n.avg))
  dat$B     <- rep(1:K, each = n_jk, times = J)

  # 每层的处理效应（使用重尾分布产生异质性）
  tau_a0_k <- rt(K, df = 3)
  tau_a1_k <- rt(K, df = 3) + tau_a0_k

  # 潜在结局 (mechanism a=0 / a=1)
  dat$y_d0_a0 <- rnorm(N)
  dat$y_d1_a0 <- dat$y_d0_a0 + tau_a0_k[dat$B] + rnorm(N)
  dat$y_d0_a1 <- rnorm(N)
  dat$y_d1_a1 <- dat$y_d0_a1 + tau_a1_k[dat$B] + rnorm(N)

  # 依从概率（按层，正数保证单调性）
  d_a0 <- abs(rt(K, df = 3)) * 0.5
  d_a1 <- abs(rt(K, df = 3)) * 0.5 + d_a0

  # 潜在中介 D (连续 -> 二值)
  dat$d_z0_a0 <- rnorm(N) * 0.1
  dat$d_z1_a0 <- dat$d_z0_a0 + d_a0[dat$B] + abs(rnorm(N) * 0.1)
  dat$d_z0_a1 <- rnorm(N) * 0.1
  dat$d_z1_a1 <- dat$d_z0_a1 + d_a1[dat$B] + abs(rnorm(N) * 0.1)

  dat$d_z0_a0 <- ifelse(dat$d_z0_a0 > 0, 1, 0)
  dat$d_z1_a0 <- ifelse(dat$d_z1_a0 > 0, 1, 0)
  dat$d_z0_a1 <- ifelse(dat$d_z0_a1 > 0, 1, 0)
  dat$d_z1_a1 <- ifelse(dat$d_z1_a1 > 0, 1, 0)

  # 观测结局（由 Z 和 A 决定的潜在结局选择）
  dat$y_z1_a0 <- ifelse(dat$d_z1_a0 == 1, dat$y_d1_a0, dat$y_d0_a0)
  dat$y_z0_a0 <- ifelse(dat$d_z0_a0 == 1, dat$y_d1_a0, dat$y_d0_a0)
  dat$y_z1_a1 <- ifelse(dat$d_z1_a1 == 1, dat$y_d1_a1, dat$y_d0_a1)
  dat$y_z0_a1 <- ifelse(dat$d_z0_a1 == 1, dat$y_d1_a1, dat$y_d0_a1)

  dat <- as.data.frame(dat)

  # 计算真实 CADE（仅 complier 的平均直接效应）
  cdata_a0 <- dat[dat$d_z1_a0 - dat$d_z0_a0 == 1, ]
  cdata_a1 <- dat[dat$d_z1_a1 - dat$d_z0_a1 == 1, ]
  cade_a0  <- mean(cdata_a0$y_z1_a0 - cdata_a0$y_z0_a0)
  cade_a1  <- mean(cdata_a1$y_z1_a1 - cdata_a1$y_z0_a1)

  cat("True CADE(0):", cade_a0, "\n")
  cat("True CADE(1):", cade_a1, "\n")

  list(data = dat, N = N, n.avg = n.avg, cade_a0 = cade_a0, cade_a1 = cade_a1)
}

# gen_data.R
# 依赖全局变量: J, K, n_jk
# 输出全局变量: data, N, n.avg, cade_a0, cade_a1

n.avg <- n_jk * K
N <- J * n.avg

data <- list()
data$id <- as.factor(rep(1:J, each = n.avg))
data$B <- rep(1:K, each = n_jk, times = J)

tau_a0_k <- rt(K, df = 3)
tau_a1_k <- rt(K, df = 3) + tau_a0_k

data$y_d0_a0 <- rnorm(N)
data$y_d1_a0 <- data$y_d0_a0 + tau_a0_k[data$B] + rnorm(N)
data$y_d0_a1 <- rnorm(N)
data$y_d1_a1 <- data$y_d0_a1 + tau_a1_k[data$B] + rnorm(N)

d_a0 <- abs(rt(K, df = 3)) * 0.5
d_a1 <- abs(rt(K, df = 3)) * 0.5 + d_a0
data$d_z0_a0 <- rnorm(N) * 0.1
data$d_z1_a0 <- data$d_z0_a0 + d_a0[data$B] + abs(rnorm(N) * 0.1)
data$d_z0_a1 <- rnorm(N) * 0.1
data$d_z1_a1 <- data$d_z0_a1 + d_a1[data$B] + abs(rnorm(N) * 0.1)

data$d_z0_a0 <- ifelse(data$d_z0_a0 > 0, 1, 0)
data$d_z1_a0 <- ifelse(data$d_z1_a0 > 0, 1, 0)
data$d_z0_a1 <- ifelse(data$d_z0_a1 > 0, 1, 0)
data$d_z1_a1 <- ifelse(data$d_z1_a1 > 0, 1, 0)

data$y_z1_a0 <- ifelse(data$d_z1_a0 == 1, data$y_d1_a0, data$y_d0_a0)
data$y_z0_a0 <- ifelse(data$d_z0_a0 == 1, data$y_d1_a0, data$y_d0_a0)
data$y_z1_a1 <- ifelse(data$d_z1_a1 == 1, data$y_d1_a1, data$y_d0_a1)
data$y_z0_a1 <- ifelse(data$d_z0_a1 == 1, data$y_d1_a1, data$y_d0_a1)

data <- as.data.frame(data)

# 计算真实 CADE
cdata_a0 <- data[data$d_z1_a0 - data$d_z0_a0 == 1, ]
cdata_a1 <- data[data$d_z1_a1 - data$d_z0_a1 == 1, ]
cade_a0 <- mean(cdata_a0$y_z1_a0 - cdata_a0$y_z0_a0)
cade_a1 <- mean(cdata_a1$y_z1_a1 - cdata_a1$y_z0_a1)

cat("True CADE(0):", cade_a0, "\n")
cat("True CADE(1):", cade_a1, "\n")

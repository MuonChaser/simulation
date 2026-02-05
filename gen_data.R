data <- list()

K=4
N=J * n.avg
data$id <- as.factor(rep(1:J, each = n.avg))
data$B <- rep(1:K, each = n.avg/K, times = J)
# 按0.1, 0.2, 0.3, 0.4比例分层
# strata_props <- c(0.1, 0.2, 0.3, 0.4)
# strata_sizes <- round(n.avg * strata_props)
# data$B <- rep(rep(1:K, times = strata_sizes), times = J)



tau_a0_k <- rt(K, df = 3)
tau_a1_k <- rt(K, df = 3) + tau_a0_k

data$y_d0_a0 <- rnorm(N)
data$y_d1_a0 <- data$y_d0_a0 + tau_a0_k[data$B] + rnorm(N)
data$y_d0_a1 <- rnorm(N)
data$y_d1_a1 <- data$y_d0_a1 + tau_a1_k[data$B] + rnorm(N)

# 没想好起什么名字
d_a0  <- abs(rt(K, df = 3)) * 0.05
d_a1  <- abs(rt(K, df = 3)) * 0.05 + d_a0
# d_a1  <- d_a0
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



data = as.data.frame(data)
cdata_a0 = data[data$d_z1_a0 - data$d_z0_a0 == 1, ]
cdata_a1 = data[data$d_z1_a1 - data$d_z0_a1 == 1, ]
cade_a0 = mean(cdata_a0$y_z1_a0 - cdata_a0$y_z0_a0)
cade_a1 = mean(cdata_a1$y_z1_a1 - cdata_a1$y_z0_a1)
library('mvtnorm')
library('tidyverse')

# set.seed(202501022)
DEBUG = TRUE
#### 核心辅助函数 ####

## 分类随机抽样
## 输出：1,2,2,1,1,1,3,...
rcategorical = function(n, p){
  rv <- runif(n, 0, 1)
  out <- rep(1, n)
  l <- length(p) - 1
  for(i in 1:l){
    out <- ifelse(rv > sum(p[1:i]), i+1, out)
  }
  return(out)
}



## 生成类大小（统一函数）
generate_cluster_sizes = function(J, n.avg, type = "equal"){
  if(type == "equal"){
    n <- rep(n.avg, J)
  } else if(type == "unequal"){
    n  <- rep(n.avg, J) * 1:J
  }
  return(n)
}


## 生成层级分配（统一函数）
generate_strata = function(n, K=5, type = "equal"){
  N <- sum(n)
  J <- length(n)

  # 确定概率
  if(type == "equal"){
    prob <- rep(1/K, K)
  } else if(type == "proportional"){
    prob <- (K:1) / sum(K:1)
  }

  # 分配层级
  B <- rep(0, N)
  cum.n <- c(0, cumsum(n))

  for(j in 1:J){
    B[(cum.n[j]+1):cum.n[j+1]] <- sample(1:K, n[j], replace = TRUE, prob = prob)
  }

  return(B)
}


#### 核心模拟函数 ####

## 1. 生成基础数据结构
create_dataset = function(J, n.avg, K,
                       cluster_type = "equal",
                       strata_type = "equal"){
  n <- generate_cluster_sizes(J, n.avg, cluster_type)
  B <- generate_strata(n, K, strata_type)
  id <- rep(1:J, n)

  data <- data.frame(id = id, B = B)

  return(list(data = data, n = n, J = J, K = K))
}


## 2. 生成所有潜在结果
generate_potential_outcomes = function(sim_data,
                                     # D 的参数
                                     d_mean_base = 0.4, d_sd_base = 0.15,
                                     d_mean_delta_z = 0.3, d_sd_delta_z = 0.08,
                                     d_mean_delta_a = 0.1, d_sd_delta_a = 0.05,
                                     d_threshold = 0.5,
                                     # Y 的参数
                                     y_mean_base = 10, y_sd_base = 2,
                                     y_mean_delta_d = 5, y_sd_delta_d = 1,
                                     y_mean_delta_a = 2, y_sd_delta_a = 1,
                                     interaction_mean = 2, interaction_sd = 0.5,
                                     strata_effect = c(1,2,3,4,5)
                                     ){

  data <- sim_data$data
  N <- nrow(data)

  # --- A. 为 D (依从行为) 生成潜在结果 ---
  # D是关于z和a的函数，D = f(z, a)

  # a) 生成基础得分和变化量
  d_base_score <- rnorm(N, d_mean_base, d_sd_base)
  d_delta_z <- pmax(0.01, rnorm(N, d_mean_delta_z, d_sd_delta_z)) # 确保z=1的效果为正
  d_delta_a <- rnorm(N, d_mean_delta_a, d_sd_delta_a)

  # b) 基于得分和阈值，直接计算4种组合下的潜在D值(0/1)
  data$d_z0_a0 <- as.integer(d_base_score > d_threshold)
  data$d_z0_a1 <- as.integer((d_base_score + d_delta_a) > d_threshold)
  data$d_z1_a0 <- as.integer((d_base_score + d_delta_z) > d_threshold)
  data$d_z1_a1 <- as.integer((d_base_score + d_delta_z + d_delta_a) > d_threshold)

  # --- B. 为 Y (结果) 生成潜在结果 ---
  # Y是关于d和a的函数, Y = f(d, a)

  # a) 生成基础结果和变化量
  y_base <- rnorm(N, y_mean_base, y_sd_base) # 这就是 Y(d=0, a=0)
  

  y_delta_d <- rnorm(N, y_mean_delta_d, y_sd_delta_d) + data$B
  y_delta_a <- rnorm(N, y_mean_delta_a, y_sd_delta_a)

  # b) 计算4种组合下的潜在结果
  data$y_d0_a0 <- y_base
  data$y_d0_a1 <- y_base + y_delta_a
  data$y_d1_a0 <- y_base + y_delta_d
  data$y_d1_a1 <- y_base + y_delta_d + y_delta_a + rnorm(N, interaction_mean, interaction_sd)

  sim_data$data <- data
  return(sim_data)
}


#### 分析函数 ####

## 计算CADE (Complier Average Direct Effect)
calculate_cade = function(sim_data) {
  data <- sim_data$data
  ded_a0 <- mean(data$d_z1_a0 - data$d_z0_a0)
  ded_a1 <- mean(data$d_z1_a1 - data$d_z0_a1)
  dey_a1 <- mean(data$y_d1_a1 - data$y_d0_a1)
  dey_a0 <- mean(data$y_d1_a0 - data$y_d0_a0)

  sim_data$real$ded_a0 <- ded_a0
  sim_data$real$ded_a1 <- ded_a1
  sim_data$real$dey_a0 <- dey_a0
  sim_data$real$dey_a1 <- dey_a1


  # --- CADE for a=0 ---
  is_complier_a0 <- data$d_z1_a0 == 1 & data$d_z0_a0 == 0
  effects_a0 <- data$y_d1_a0[is_complier_a0] - data$y_d0_a0[is_complier_a0]
  cade_a0 <- mean(effects_a0, na.rm = TRUE)

  # --- CADE for a=1 ---
  is_complier_a1 <- data$d_z1_a1 == 1 & data$d_z0_a1 == 0
  effects_a1 <- data$y_d1_a1[is_complier_a1] - data$y_d0_a1[is_complier_a1]
  cade_a1 <- mean(effects_a1, na.rm = TRUE)

  # 将结果附加到sim_data对象并返回
  sim_data$real$CADE_a0 <- cade_a0
  sim_data$real$CADE_a1 <- cade_a1

  if (DEBUG) {
    cat("Count of compliers (a=0):", sum(is_complier_a0), "\n")
    cat("Count of compliers (a=1):", sum(is_complier_a1), "\n")
    cat("Count of always-takers (a=0):", sum(data$d_z1_a0 == 1 & data$d_z0_a0 == 1), "\n")
    cat("Count of always-takers (a=1):", sum(data$d_z1_a1 == 1 & data$d_z0_a1 == 1), "\n")
    cat("Count of never-takers (a=0):", sum(data$d_z1_a0 == 0 & data$d_z0_a0 == 0), "\n")
    cat("Count of never-takers (a=1):", sum(data$d_z1_a1 == 0 & data$d_z0_a1 == 0), "\n")
  }

  return(sim_data)
}


## 3. 分配处理和类别
assign_cluster_level = function(sim_data, prop_high = 0.5){
  J <- sim_data$J
  n_high <- round(J * prop_high)
  a <- rep(0, J)
  high_idx <- sample(1:J, n_high, replace = FALSE)
  a[high_idx] <- 1

  # 扩展为个体水平
  sim_data$data$a <- rep(a, sim_data$n)
  sim_data$a <- a
  return(sim_data)
}


## 4. 分配个体处理
assign_treatment_group = function(sim_data, prop_treat_high = 0.5, prop_treat_low = 0.5){
  data <- sim_data$data
  n <- sim_data$n
  J <- sim_data$J
  K <- sim_data$K
  N <- nrow(data)
  a <- sim_data$a

  z <- rep(0, N)
  cum.n <- c(0, cumsum(n))

  for(j in 1:J){
    idx_start <- cum.n[j] + 1
    idx_end <- cum.n[j + 1]
    idx_class <- idx_start:idx_end

    prop_treat <- ifelse(a[j] == 1, prop_treat_high, prop_treat_low)

    # 在每个层内精确分配
    for(k in 1:K){
      idx_stratum <- idx_class[data$B[idx_class] == k]
      n_stratum <- length(idx_stratum)

      if(n_stratum > 0){
        n_treat <- round(n_stratum * prop_treat)
        if(n_treat > 0){
          treat_idx <- sample(idx_stratum, n_treat, replace = FALSE)
          z[treat_idx] <- 1
        }
      }
    }
  }

  sim_data$data$z <- z
  return(sim_data)
}


## 5. 揭示观测结果
cal_outcomes = function(sim_data){
  data <- sim_data$data

  # --- A. 揭示依从行为 D ---
  # 根据分配的 a 和 z，选择对应的潜在D值
  data$D <- case_when(
    data$a == 0 & data$z == 0 ~ data$d_z0_a0,
    data$a == 1 & data$z == 0 ~ data$d_z0_a1,
    data$a == 0 & data$z == 1 ~ data$d_z1_a0,
    data$a == 1 & data$z == 1 ~ data$d_z1_a1
  )

  # --- B. 揭示结果 Y ---
  # 根据实现的 D 和分配的 a，选择对应的潜在结果
  data$y <- case_when(
    data$D == 0 & data$a == 0 ~ data$y_d0_a0,
    data$D == 1 & data$a == 0 ~ data$y_d1_a0,
    data$D == 0 & data$a == 1 ~ data$y_d0_a1,
    data$D == 1 & data$a == 1 ~ data$y_d1_a1
  )

  sim_data$data <- data
  return(sim_data)
}

## 估计 DED (Direct Effect of D)
## 使用加权最小二乘法，分别估计 a=0 和 a=1 组的 CADE
estimate_ded = function(sim_data){
  data <- sim_data$data
  N <- nrow(data)      # 总样本量
  K <- sim_data$K      # 层数
  J <- sim_data$J      # Cluster总数

  for(k in 2:K){
    n_k <- sum(data$B == k)  # 第k层的人数
    data[[paste0("pi_star_", k)]] <- (n_k / N) * as.integer(data$B == k)
  }

  # --------------------------------------------------------------
  # 步骤 2: 计算所需的统计量
  # --------------------------------------------------------------
  # n_k: 每层总人数
  n_k_table <- data %>% count(B, name = "n_k")

  # n_jk: 每个cluster在每层的人数
  cluster_strata_stats <- data %>%
    group_by(id, a, B) %>%
    summarise(n_jk = n(), .groups = "drop")

  # n_jkz
  cluster_strata_z_stats <- data %>%
    group_by(id, a, B, z) %>%
    summarise(n_jkz = n(), .groups = "drop")

  # J_a: 每个a组的cluster数
  J_a_table <- data %>%
    distinct(id, a) %>%
    count(a, name = "J_a")


  data <- data %>%
    left_join(n_k_table, by = "B") %>%
    left_join(cluster_strata_stats, by = c("id", "a", "B")) %>%
    left_join(J_a_table, by = "a") %>%
    left_join(cluster_strata_z_stats, by = c("id", "a", "B", "z")) %>%
    mutate(
      # 转换后的因变量: D* = n_jk * J * D / n_k
      D_star = D * n_jk * J  / n_k,
      Y_star = y * n_jk * J  / n_k,
      # 回归权重: w = 1 / (J_a * n_jkz)
      # n_jkz 是该个体所在 (cluster j, 层 k, z值) 组合的人数
      weight = 1 / (J_a * n_jkz)
    )

  # --------------------------------------------------------------
  # 步骤 4: 按 a 分组
  # --------------------------------------------------------------
  data_a0 <- data %>% filter(a == 0)  # 低资源组
  data_a1 <- data %>% filter(a == 1)  # 高资源组

 
  pi_vars <- paste0("pi_star_", 2:K)  # 从第2层开始，第1层作为基准

  # 构建公式：主效应 + 交互效应
  formula_str_d <- paste0(
    "D_star ~ z + ",                               # z 的主效应
    paste(pi_vars, collapse = " + "),              # 层级的主效应
    " + ",
    paste0("z:", pi_vars, collapse = " + ")        # z 与各层的交互效应
  )

  # 在 a=0 组运行 WLS 回归
  model_d_a0 <- lm(as.formula(formula_str_d), data = data_a0, weights = weight)
  ded_a0 <- coef(model_d_a0)["z"]
  
  # 在 a=1 组运行 WLS 回归
  model_d_a1 <- lm(as.formula(formula_str_d), data = data_a1, weights = weight)
  ded_a1 <- coef(model_d_a1)["z"]

  formula_str_y <- paste0(
    "Y_star ~ z + ",                               # z 的主效应
    paste(pi_vars, collapse = " + "),              # 层级的主效应
    " + ",
    paste0("z:", pi_vars, collapse = " + ")        # z 与各层的交互效应
  )

  model_y_a0 <- lm(as.formula(formula_str_y), data = data_a0, weights = weight)
  dey_a0 <- coef(model_y_a0)["z"]

  model_y_a1 <- lm(as.formula(formula_str_y), data = data_a1, weights = weight)
  dey_a1 <- coef(model_y_a1)["z"]



  sim_data$simu$data <- data
  sim_data$simu$ded_a0 <- ded_a0  
  sim_data$simu$ded_a1 <- ded_a1  
  sim_data$simu$dey_a0 <- dey_a0  
  sim_data$simu$dey_a1 <- dey_a1  

  # 调试模式下保存详细信息
  if(DEBUG){
    sim_data$simu$data_a0 <- data_a0
    sim_data$simu$data_a1 <- data_a1
    sim_data$simu$model_d_a0 <- model_d_a0
    sim_data$simu$model_d_a1 <- model_d_a1
    sim_data$simu$model_y_a0 <- model_y_a0
    sim_data$simu$model_y_a1 <- model_y_a1
  }

  return(sim_data)
}


#### 抽象测试函数 ####

#' 研究无偏性
#' 
#' @param sim_data_base 基础模拟数据对象（包含潜在结果）
#' @param pipeline_fun 模拟流程函数（接受sim_data返回sim_data）
#' @param n_sims 模拟次数
#' @param get_est 从sim_data中提取估计量的函数
#' @param get_true 从sim_data中提取真实值的函数
#' @return 包含估计值、偏差等信息的列表
run_unbiasedness_test <- function(sim_data_base, pipeline_fun, n_sims = 200, 
                                  get_est = function(x) x$simu$ded_a0, 
                                  get_true = function(x) x$real$ded_a0) {
  
  # 获取基准真实值（假设针对该有限人口是固定的）
  true_val <- get_true(sim_data_base)
  
  cat(sprintf("\n--- 开始无偏性测试 (N_sims = %d) ---\n", n_sims))
  
  estimates <- replicate(n_sims, {
    sim_res <- pipeline_fun(sim_data_base)
    get_est(sim_res)
  })
  
  mean_est <- mean(estimates)
  bias <- mean_est - true_val
  mse <- mean((estimates - true_val)^2)
  
  cat(sprintf("真实值 (Truth): %.4f\n", true_val))
  cat(sprintf("估计均值 (Mean Est.): %.4f\n", mean_est))
  cat(sprintf("偏差 (Bias): %.4f\n", bias))
  cat(sprintf("均方误差 (MSE): %.4f\n", mse))
  
  # 简单绘图
  hist(estimates, breaks = 20, main = "估计量分布", xlab = "Estimate Value", col = "lightblue", border = "white")
  abline(v = true_val, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c("Estimates", "Truth"), fill = c("lightblue", "red"), density = c(NA, NA), border = c("white", "white"), lty = c(NA, 2), lwd = c(NA, 2))
  
  return(list(
    estimates = estimates,
    truth = true_val,
    bias = bias,
    mse = mse
  ))
}


#' 研究渐进特性 (Consistency/Asymptotic)
#'
#' @param J_values 样本量向量（Cluster数量）
#' @param setup_factory_fun 接受参数 J 并返回 sim_data_base 的函数
#' @param pipeline_fun 模拟流程函数
#' @param n_sims_per_j 每个样本量下的模拟次数
#' @param get_est 提取估计量函数
#' @param get_true 提取真实值函数
#' @return 包含所有模拟结果的数据框
run_asymptotic_test <- function(J_values, setup_factory_fun, pipeline_fun, n_sims_per_j = 100,
                                get_est = function(x) x$simu$ded_a0,
                                get_true = function(x) x$real$ded_a0) {
  
  all_results <- data.frame()
  
  cat("\n--- 开始渐进特性测试 ---\n")
  
  for(J in J_values) {
    cat(sprintf("Processing J = %d...\n", J))
    
    # 1. 生成该样本量下的总体
    sim_data_base <- setup_factory_fun(J)
    true_val <- get_true(sim_data_base)
    
    # 2. 运行多次模拟
    estimates <- replicate(n_sims_per_j, {
      sim_res <- pipeline_fun(sim_data_base)
      get_est(sim_res)
    })
    
    # 3. 收集结果
    tmp_df <- data.frame(
      J = J,
      estimate = estimates,
      truth = true_val,
      bias = estimates - true_val
    )
    all_results <- rbind(all_results, tmp_df)
  }
  
  # 4. 汇总与绘图
  summary_stats <- all_results %>%
    group_by(J) %>%
    summarise(
      mean_bias = mean(bias),
      mse = mean(bias^2),
      sd_bias = sd(bias)
    )
  print(summary_stats)
  
  p <- ggplot(all_results, aes(x = factor(J), y = bias)) +
    geom_boxplot(fill = "lightgreen", alpha = 0.7) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
    labs(
      title = "估计量偏差的渐进特性",
      subtitle = "Bias = Estimate - Truth",
      x = "Cluster 数量 (J)",
      y = "偏差 (Bias)"
    ) +
    theme_minimal()
  
  print(p)
  
  return(all_results)
}



# 设定 DEBUG = FALSE 以在大量模拟中提高效率
DEBUG = FALSE

# --- 定义通用的 pipeline ---
my_pipeline <- function(data) {
  data %>%
    assign_cluster_level(prop_high = 0.5) %>%
    assign_treatment_group(prop_treat_high = 0.7, prop_treat_low = 0.3) %>%
    cal_outcomes() %>%
    estimate_ded()
}

# --- 示例 1: 单一样本量的无偏性测试 ---
cat("\n========== 示例 1: 研究无偏性 (J=50) ==========\n")

# 创建基础数据
base_data <- create_dataset(J = 30, n.avg = 20, K = 3) %>%
  generate_potential_outcomes() %>%
  calculate_cade()

# 运行测试：研究 ded_a0
res_unbiased <- run_unbiasedness_test(
  sim_data_base = base_data,
  pipeline_fun = my_pipeline,
  n_sims = 2000,
  get_est = function(x) x$simu$ded_a0,
  get_true = function(x) x$real$ded_a0
)


# # --- 示例 2: 渐进特性测试 (不同J值) ---
# cat("\n========== 示例 2: 研究渐进特性 ==========\n")

# # 定义如何生成不同J的数据
# my_setup_factory <- function(j) {
#   create_dataset(J = j, n.avg = 60, K = 3) %>%
#     generate_potential_outcomes() %>%
#     calculate_cade()
# }

# # 运行测试
# res_asymp <- run_asymptotic_test(
#   J_values = c(20, 50, 100, 200),
#   setup_factory_fun = my_setup_factory,
#   pipeline_fun = my_pipeline,
#   n_sims_per_j = 100,
#   get_est = function(x) x$simu$ded_a0,
#   get_true = function(x) x$real$ded_a0
# )

# # 保存图片
# ggsave("simulation_asymptotic_bias.png", width = 8, height = 6)




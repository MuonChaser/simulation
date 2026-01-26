library('mvtnorm')
library('tidyverse')

set.seed(202501022)

#### 生成基础数据结构
#### J个类，每个类n_j个个体，每个个体属于层B_ij (1到K)


## 分类随机抽样函数
rcategorical = function(n, p){
  rv <- runif(n, 0, 1)
  out <- rep(1, n)
  l <- length(p) - 1
  for(i in 1:l){
    out <- ifelse(rv > sum(p[1:i]), i+1, out)
  }
  return(out)
}


GenClusterSize_equal = function(J, n.avg){
  n <- rep(n.avg, J)
  return(n)
}


GenClusterSize_unequal = function(J, n.avg){
  n.index <- rcategorical(J, c(0.25, 0.1, 0.1, 0.1, 0.2, 0.25))
  n <- round(n.avg * c(0.5, 0.75, 1, 1.5, 2, 2.5)[n.index])
  n <- pmax(n, 1)
  return(n)
}


## 生成类大小 - 极端分布
GenClusterSize_extreme = function(J, n.avg){
  n.index <- rcategorical(J, c(0.7, 0.2, 0.1))
  n <- round(n.avg * c(0.3, 1, 5)[n.index])
  n <- pmax(n, 1)
  return(n)
}




## 生成层级分配 - 等概率分配（均匀分布）
GenStrata_equal = function(n, K){
  N <- sum(n)
  J <- length(n)
  prob <- rep(1/K, K)
  
  B <- rep(0, N)
  cum.n <- c(0, cumsum(n))
  
  for(j in 1:J){
    B[(cum.n[j]+1):cum.n[j+1]] <- sample(1:K, n[j], replace = TRUE, prob = prob)
  }
  
  return(B)
}


## 生成层级分配 - 等比例分配（递减）
GenStrata_proportional = function(n, K){
  N <- sum(n)
  J <- length(n)
  prob <- (K:1)
  prob <- prob/sum(prob)
  
  B <- rep(0, N)
  cum.n <- c(0, cumsum(n))
  
  for(j in 1:J){
    B[(cum.n[j]+1):cum.n[j+1]] <- sample(1:K, n[j], replace = TRUE, prob = prob)
  }
  
  return(B)
}


## 生成完整数据结构
GenBaseData = function(n, K, gen_strata_func, ...){
  # n: 每个类的大小向量
  # K: 层数
  # gen_strata_func: 生成层级分配的函数
  # ...: 传递给gen_strata_func的额外参数
  
  J <- length(n)
  
  # 生成层级分配
  B <- gen_strata_func(n, K, ...)
  
  # 生成类ID
  id <- rep(1:J, n)
  
  # 组合成数据框
  data <- data.frame(
    id = id,
    B = B
  )
  
  return(list(data = data, n = n, J = J, K = K))
}


#### 分配处理向量 ####

## 第一步：给每个类分配高/低水平 (a_j = 0 或 1)
## prop_high: 高水平组的比例，精确分配
AssignClassLevel = function(J, prop_high = 0.5){
  n_high <- round(J * prop_high)  # 高水平类的数量
  a <- rep(0, J)
  high_idx <- sample(1:J, n_high, replace = FALSE)  # 随机选择哪些类是高水平
  a[high_idx] <- 1
  return(a)
}


## 第二步：在每个类的每个层内分配处理 (z = 0 或 1)
## base_data: GenBaseData返回的list
## a: 类水平向量 (长度为J)
## prop_treat_high: 高水平组(a=1)中的处理比例
## prop_treat_low: 低水平组(a=0)中的处理比例
AssignTreatment = function(base_data, a, prop_treat_high = 0.5, prop_treat_low = 0.5){
  data <- base_data$data
  n <- base_data$n
  J <- base_data$J
  K <- base_data$K
  N <- nrow(data)

  z <- rep(0, N)
  cum.n <- c(0, cumsum(n))

  for(j in 1:J){
    # 获取该类的索引范围
    idx_start <- cum.n[j] + 1
    idx_end <- cum.n[j + 1]
    idx_class <- idx_start:idx_end

    # 确定该类的处理比例
    prop_treat <- ifelse(a[j] == 1, prop_treat_high, prop_treat_low)

    # 在每个层内精确分配
    for(k in 1:K){
      # 找到该类该层的所有个体索引
      idx_stratum <- idx_class[data$B[idx_class] == k]
      n_stratum <- length(idx_stratum)

      if(n_stratum > 0){
        # 精确计算处理组人数
        n_treat <- round(n_stratum * prop_treat)
        if(n_treat > 0){
          treat_idx <- sample(idx_stratum, n_treat, replace = FALSE)
          z[treat_idx] <- 1
        }
      }
    }
  }

  return(z)
}


## 完整的处理分配函数（包装前两步）
AssignAllTreatment = function(base_data, prop_high = 0.5,
                               prop_treat_high = 0.5, prop_treat_low = 0.5){
  J <- base_data$J

  # 第一步：分配类水平
  a <- AssignClassLevel(J, prop_high)

  # 第二步：分配个体处理
  z <- AssignTreatment(base_data, a, prop_treat_high, prop_treat_low)

  # 将a扩展为个体水平的向量（方便后续使用）
  a_individual <- rep(a, base_data$n)

  # 添加到数据中
  base_data$data$a <- a_individual
  base_data$data$z <- z
  base_data$a <- a  # 类水平的a向量

  return(base_data)
}

#### 模拟不遵从行为（Non-compliance）####

## 生成个体初始处理倾向
## 每个个体有两个基础值: d0 (被分配到对照时) 和 d1 (被分配到处理时)
## 设计原则：
##   - 使用"基础倾向 + 偏移"的方式，确保 d1 > d0
##   - base_propensity: 个体的基础处理倾向（个体异质性）
##   - gap: d1与d0之间的差距（分配效应强度）
GenBasePropensity = function(N, mean_base = 0.5, sd_base = 0.15,
                              mean_gap = 0.3, sd_gap = 0.08){
  # 生成基础倾向（个体的固有处理偏好）
  base_propensity <- rnorm(N, mean_base, sd_base)

  # 生成d0和d1之间的差距（确保为正）
  gap <- pmax(0.1, rnorm(N, mean_gap, sd_gap))

  # 计算d0和d1，以base为中心对称分布
  d0 <- base_propensity - gap / 2
  d1 <- base_propensity + gap / 2

  # 限制在合理范围内 [0.02, 0.98]
  d0 <- pmax(0.02, pmin(0.98, d0))
  d1 <- pmax(0.02, pmin(0.98, d1))

  # 确保 d1 > d0（处理边界情况）
  swap_idx <- d1 <= d0
  if(any(swap_idx)){
    temp <- d1[swap_idx]
    d1[swap_idx] <- d0[swap_idx] + 0.1
    d1 <- pmin(0.98, d1)
  }

  return(data.frame(d0 = d0, d1 = d1))
}


## 计算实际处理状态（含同伴效应的迭代）
## 设计原则：
##   - 使用上一轮的D值计算同伴效应（同步更新，避免顺序依赖）

SimulateCompliance = function(sim_data,
                               # 基础倾向参数
                               mean_base = 0.5, sd_base = 0.15,
                               mean_gap = 0.3, sd_gap = 0.08,
                               # 同伴效应参数
                               gamma = 0.2,
                               # 迭代参数
                               threshold = 0.5,
                               n_iter = 2,
                               verbose = TRUE){

  data <- sim_data$data
  n <- sim_data$n
  J <- sim_data$J
  N <- nrow(data)
  cum.n <- c(0, cumsum(n))

  # 生成基础倾向值
  propensity <- GenBasePropensity(N, mean_base, sd_base, mean_gap, sd_gap)
  data$d0 <- propensity$d0
  data$d1 <- propensity$d1

  # 根据分配确定初始处理值（无同伴效应）
  data$treat_value <- ifelse(data$z == 1, data$d1, data$d0)

  # 初始实际处理状态
  data$D <- as.integer(data$treat_value > threshold)

  # 预计算每个个体所属的类
  class_idx_list <- vector("list", J)
  for(j in 1:J){
    class_idx_list[[j]] <- (cum.n[j] + 1):cum.n[j + 1]
  }

  # 迭代更新（同伴效应）- 同步更新固定次数
  for(iter in 1:n_iter){
    D_old <- data$D

    # 新的处理值（基于上一轮D计算）
    treat_value_new <- rep(0, N)

    for(j in 1:J){
      idx_class <- class_idx_list[[j]]
      n_class <- length(idx_class)

      # 该类中上一轮接受处理的比例
      class_D <- D_old[idx_class]

      for(i in idx_class){
        # 计算同伴比例（排除自己）
        i_local <- which(idx_class == i)
        others_D <- class_D[-i_local]
        peer_prop <- mean(others_D)

        # 基础值
        base_val <- ifelse(data$z[i] == 1, data$d1[i], data$d0[i])

        # 同伴效应：偏离0.5的部分影响处理值
        peer_effect <- gamma * (peer_prop - 0.5)

        # 更新处理值
        treat_value_new[i] <- pmax(0, pmin(1, base_val + peer_effect))
      }
    }

    # 同步更新所有处理值
    data$treat_value <- treat_value_new

    # 更新实际处理状态
    data$D <- as.integer(data$treat_value > threshold)

    if(verbose) cat("完成第", iter, "轮迭代\n")
  }

  sim_data$data <- data
  sim_data$n_iter <- n_iter
  return(sim_data)
}


#### 生成结果变量 (y0, y1, y) ####
GenerateOutcomes = function(sim_data,
                            base_outcome = 10,
                            strata_effects = c(1, 2, 5), # K=3 示例
                            error_sd = 2){

  data <- sim_data$data
  K <- sim_data$K
  N <- nrow(data)

  # 确保 strata_effects 长度与 K 匹配
  if(length(strata_effects) != K){
    stop("strata_effects 的长度必须等于 K")
  }

  # 1. 生成 y0
  # y0 受层级 B 和个体随机误差影响
  data$y0 <- base_outcome + data$B * 0.5 + rnorm(N, 0, error_sd)

  # 2. 生成 y1，确保 y1 >= y0
  # y1 - y0 的差值由层决定
  delta_k <- strata_effects[data$B]
  data$y1 <- data$y0 + delta_k

  # 3. 生成观察到的 y
  # 根据用户要求：d=1 时 y=y0，否则 y=y1 (d用D表示)
  data$y <- ifelse(data$D == 1, data$y0, data$y1)

  sim_data$data <- data
  return(sim_data)
}


## 完整的模拟函数（从分配到不遵从）
RunFullSimulation = function(base_data,
                              # 分配参数
                              prop_high = 0.5,
                              prop_treat_high = 0.7,
                              prop_treat_low = 0.3,
                              # 不遵从参数 - 基础倾向
                              mean_base = 0.5, sd_base = 0.15,
                              mean_gap = 0.3, sd_gap = 0.08,
                              # 不遵从参数 - 同伴效应
                              gamma = 0.2,
                              threshold = 0.5,
                              n_iter = 2,
                              # 结果参数
                              base_outcome = 10,
                              strata_effects = c(1, 2, 5),
                              error_sd = 2,
                              verbose = TRUE){

  # 分配处理
  sim_data <- AssignAllTreatment(base_data, prop_high, prop_treat_high, prop_treat_low)

  # 模拟不遵从
  sim_data <- SimulateCompliance(sim_data,
                                  mean_base, sd_base, mean_gap, sd_gap,
                                  gamma, threshold, n_iter, verbose)
                                  
  # 生成结果
  sim_data <- GenerateOutcomes(sim_data, base_outcome, strata_effects, error_sd)

  return(sim_data)
}


#### 测试不遵从模拟 ####

cat("========================================\n")
cat("       完整模拟测试 (含结果)\n")
cat("========================================\n\n")

n <- GenClusterSize_equal(20, 30)
# base_data <- GenBaseData(n, K = 3, gen_strata_func = GenStrata_equal)
base_data <- GenBaseData(n, K = 3, gen_strata_func = GenStrata_proportional)


# 定义结果生成参数
strata_effects_test <- c(2, 4, 6) # y1-y0 for B=1, 2, 3

sim_data <- RunFullSimulation(base_data,
                               # 分配参数
                               prop_high = 0.5,
                               prop_treat_high = 0.7,
                               prop_treat_low = 0.3,
                               # 基础倾向参数
                               mean_base = 0.4, sd_base = 0.15,
                               mean_gap = 0.3, sd_gap = 0.08,
                               # 同伴效应
                               n_iter = 2,
                               gamma = 0.2,
                               # 结果参数
                               strata_effects = strata_effects_test)

cat("\n========== 基本统计 ==========\n")
cat("总样本量:", nrow(sim_data$data), "\n")
cat("类数量:", sim_data$J, "\n")
cat("每类人数:", unique(sim_data$n)[1], "\n")

cat("\n========== 分配与实际接受 ==========\n")
cat("分配到处理组(z=1)比例:", round(mean(sim_data$data$z), 3), "\n")
cat("实际接受处理(D=1)比例:", round(mean(sim_data$data$D), 3), "\n")

# 遵从率分析
cat("\n========== 遵从率分析 ==========\n")
d <- sim_data$data
cat("总体遵从率:", round(mean(d$z == d$D), 3), "\n")
cat("  - Compliers (z=1 & D=1):", sum(d$z == 1 & d$D == 1), "\n")
cat("  - Never-takers (z=1 & D=0):", sum(d$z == 1 & d$D == 0), "\n")
cat("  - Always-takers (z=0 & D=1):", sum(d$z == 0 & d$D == 1), "\n")
cat("  - Compliers (z=0 & D=0):", sum(d$z == 0 & d$D == 0), "\n")

# 按分配组分析
cat("\n========== 按分配组分析 ==========\n")
cat("被分配处理(z=1)中实际接受(D=1):",
    round(mean(d$D[d$z == 1]), 3), "\n")
cat("被分配对照(z=0)中实际接受(D=1):",
    round(mean(d$D[d$z == 0]), 3), "\n")

# 按类水平分析
cat("\n========== 按类水平分析 ==========\n")
cat("高水平类(a=1)中:\n")
cat("  分配到处理比例:", round(mean(d$z[d$a == 1]), 3), "\n")
cat("  实际接受比例:", round(mean(d$D[d$a == 1]), 3), "\n")
cat("低水平类(a=0)中:\n")
cat("  分配到处理比例:", round(mean(d$z[d$a == 0]), 3), "\n")
cat("  实际接受比例:", round(mean(d$D[d$a == 0]), 3), "\n")

# 验证 y1 >= y0
cat("\n========== 结果变量验证 ==========\n")
cat("y0 均值:", round(mean(d$y0), 2), "\n")
cat("y1 均值:", round(mean(d$y1), 2), "\n")
cat("y 均值:", round(mean(d$y), 2), "\n")
cat("y1 >= y0 验证:", all(d$y1 >= d$y0), "\n")
cat("\n按层级的 (y1-y0) 均值:\n")
y_diff_summary <- d %>%
  group_by(B) %>%
  summarise(mean_y1_minus_y0 = mean(y1 - y0))
print(y_diff_summary)
cat("预设的 (y1-y0) 值:", strata_effects_test, "\n")


# 数据预览
cat("\n========== 数据预览 ==========\n")
print(head(d[, c("id", "B", "a", "z", "D", "y0", "y1", "y")], 15))

# 各类的处理接受情况
cat("\n========== 各类处理接受情况（前10类）==========\n")
class_summary <- d %>%
  group_by(id) %>%
  summarise(
    n = n(),
    a = first(a),
    z_prop = round(mean(z), 2),
    D_prop = round(mean(D), 2),
    .groups = "drop"
  ) %>%
  head(10)
print(class_summary)

# 保存数据
write_csv(sim_data$data, "simulated_data.csv")
cat("\n数据已保存到 simulated_data.csv\n")

#### 计算每个类、每个层的描述性统计 (按Z分组) ####
cat("\n========================================\n")
cat("    每个类-层单元的描述性统计 (按Z分组)\n")
cat("========================================\n\n")

# 根据用户最新澄清，按 (id, B) 分组，计算不同z水平下的D和Y的均值
cluster_stratum_summary_by_Z <- d %>%
  group_by(id, B) %>%
  summarise(
    n_total = n(),
    n_z0 = sum(z == 0),
    n_z1 = sum(z == 1),
    mean_D_z0 = round(mean(D[z == 0]), 3),
    mean_D_z1 = round(mean(D[z == 1]), 3),
    mean_Y_z0 = round(mean(y[z == 0]), 3),
    mean_Y_z1 = round(mean(y[z == 1]), 3),
    .groups = 'drop'
  )

cat("计算完成。预览前20行：\n")
print(head(cluster_stratum_summary_by_Z, 20))

# 保存这个汇总结果
write_csv(cluster_stratum_summary_by_Z, "cluster_stratum_summary_by_Z.csv")
cat("\n按Z分组的类-层级别汇总数据已保存到 cluster_stratum_summary_by_Z.csv\n")

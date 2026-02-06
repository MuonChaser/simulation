source("CADEreg_new.R")
source("assign.R")
library(experiment)
library(ggplot2)
library(foreach)
library(doParallel)
library(doRNG)

# source("gen_data.R")

N_simu = 10000
set.seed(123)

# 注册并行后端
n_cores <- 100*parallel::detectCores() - 1
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
  d_group  <- assign_no_strata(data)
  r_strata <- CADEreg_new(d_strata)
  r_group  <- CADEreg(d_group)
  c(r_strata$CADE0 - cade_a0,
    r_group$CADE0  - cade_a0,
    r_strata$CADE1 - cade_a1,
    r_group$CADE1  - cade_a1)
}

stopCluster(cl)

bias_strata_a0 <- results[, 1]
bias_group_a0  <- results[, 2]
bias_strata_a1 <- results[, 3]
bias_group_a1  <- results[, 4]

cat("Strata A=0 Bias Mean:", mean(bias_strata_a0), " SD:", sd(bias_strata_a0), " RMSE:", sqrt(mean(bias_strata_a0^2)))
cat(" CP:", mean(abs(bias_strata_a0) < 1.96 * sd(bias_strata_a0)), "\n")
cat(" CI length:", 2 * 1.96 * sd(bias_strata_a0), "\n")

cat("Group A=0 Bias Mean:", mean(bias_group_a0), " SD:", sd(bias_group_a0), " RMSE:", sqrt(mean(bias_group_a0^2)), "\n")
cat(" CP:", mean(abs(bias_group_a0) < 1.96 * sd(bias_group_a0)), "\n")
cat(" CI length:", 2 * 1.96 * sd(bias_group_a0), "\n")

cat("Strata A=1 Bias Mean:", mean(bias_strata_a1), " SD:", sd(bias_strata_a1), " RMSE:", sqrt(mean(bias_strata_a1^2)), "\n")
cat(" CP:", mean(abs(bias_strata_a1) < 1.96 * sd(bias_strata_a1)), "\n")
cat(" CI length:", 2 * 1.96 * sd(bias_strata_a1), "\n")

cat("Group A=1 Bias Mean:", mean(bias_group_a1), " SD:", sd(bias_group_a1), " RMSE:", sqrt(mean(bias_group_a1^2)), "\n")
cat(" CP:", mean(abs(bias_group_a1) < 1.96 * sd(bias_group_a1)), "\n")
cat(" CI length:", 2 * 1.96 * sd(bias_group_a1), "\n")


bias_data <- data.frame(
  Bias = c(bias_strata_a0, bias_group_a0, bias_strata_a1, bias_group_a1),
  Method = factor(rep(c("Strata", "Group"), each = N_simu)),
  A = factor(rep(c("A=0", "A=1"), each = 2 * N_simu))
)


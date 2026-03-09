# save_results.R
# 依赖全局变量: results, N_simu, label, K, n_jk, J
# 输出: {label}_results.rds, {label}_stats.csv, {label}_boxplot.png, {label}_violin.png

# 提取偏差向量
bias_strata_a0 <- results[, 1]
bias_group_a0 <- results[, 2]
bias_strata_a1 <- results[, 3]
bias_group_a1 <- results[, 4]

# ========== 保存 RDS ==========
saveRDS(results, paste0(label, "_results.rds"))
cat("Saved:", paste0(label, "_results.rds"), "\n")

# ========== 计算统计量 ==========
calc_stats <- function(bias, name_method, name_a) {
  data.frame(
    Method = name_method,
    A = name_a,
    Mean = mean(bias),
    SD = sd(bias),
    RMSE = sqrt(mean(bias^2)),
    CP = mean(abs(bias) < 1.96 * sd(bias)),
    CI_Length = 2 * 1.96 * sd(bias)
  )
}

stats_df <- rbind(
  calc_stats(bias_strata_a0, "Strata", "A=0"),
  calc_stats(bias_group_a0, "Group", "A=0"),
  calc_stats(bias_strata_a1, "Strata", "A=1"),
  calc_stats(bias_group_a1, "Group", "A=1")
)

write.csv(stats_df, paste0(label, "_stats.csv"), row.names = FALSE)
cat("Saved:", paste0(label, "_stats.csv"), "\n")
print(stats_df)

# ========== 准备绑图数据 ==========
bias_data <- data.frame(
  Bias = c(bias_strata_a0, bias_group_a0, bias_strata_a1, bias_group_a1),
  Method = factor(rep(c("Strata", "Group", "Strata", "Group"), each = N_simu)),
  A = factor(rep(c("A=0", "A=0", "A=1", "A=1"), each = N_simu))
)

# ========== 箱线图 ==========
g_box <- ggplot(bias_data, aes(x = Method, y = Bias, fill = Method)) +
  ylim(-5, 5) +
  geom_boxplot() +
  facet_wrap(~ A) +
  theme_minimal() +
  labs(title = paste("Bias Distribution -", label),
       subtitle = paste("K =", K, ", n_jk =", n_jk, ", J =", J),
       x = "Method",
       y = "Bias") +
  theme(legend.position = "none")

ggsave(paste0(label, "_boxplot.png"), plot = g_box, width = 8, height = 6)
cat("Saved:", paste0(label, "_boxplot.png"), "\n")

# ========== 小提琴图 ==========
g_violin <- ggplot(bias_data, aes(x = Method, y = Bias, fill = Method)) +
  ylim(-5, 5) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ A) +
  theme_minimal() +
  labs(title = paste("Bias Distribution -", label),
       subtitle = paste("K =", K, ", n_jk =", n_jk, ", J =", J),
       x = "Method",
       y = "Bias") +
  theme(legend.position = "none")

ggsave(paste0(label, "_violin.png"), plot = g_violin, width = 8, height = 6)
cat("Saved:", paste0(label, "_violin.png"), "\n")

cat("Experiment", label, "completed!\n")

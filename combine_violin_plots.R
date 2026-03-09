# 拼接old目录下的violin图
# 检查并安装所需包
if (!require(magick, quietly = TRUE)) {
  install.packages("magick", repos = "https://cloud.r-project.org")
  library(magick)
}

# 读取所有violin图
g1_K25 <- image_read("old/g1_K25_violin.png")
g1_K50 <- image_read("old/g1_K50_violin.png")
g1_K100 <- image_read("old/g1_K100_violin.png")

g2_n100 <- image_read("old/g2_n100_violin.png")
g2_n200 <- image_read("old/g2_n200_violin.png")
g2_n500 <- image_read("old/g2_n500_violin.png")

g3_J20 <- image_read("old/g3_J20_violin.png")
g3_J50 <- image_read("old/g3_J50_violin.png")
g3_J100 <- image_read("old/g3_J100_violin.png")

# 方案1: 按3x3布局拼接所有图
all_plots <- c(g1_K25, g1_K50, g1_K100,
               g2_n100, g2_n200, g2_n500,
               g3_J20, g3_J50, g3_J100)

combined_all <- image_append(image_scale(all_plots, "x400"), stack = TRUE)
combined_all <- image_append(c(
  image_append(c(g1_K25, g1_K50, g1_K100), stack = FALSE),
  image_append(c(g2_n100, g2_n200, g2_n500), stack = FALSE),
  image_append(c(g3_J20, g3_J50, g3_J100), stack = FALSE)
), stack = TRUE)

image_write(combined_all, "old/all_violin_plots_combined.png")

# 方案2: 分组拼接
# g1组 - K参数变化
g1_combined <- image_append(c(g1_K25, g1_K50, g1_K100), stack = FALSE)
image_write(g1_combined, "old/g1_violin_combined.png")

# g2组 - n参数变化
g2_combined <- image_append(c(g2_n100, g2_n200, g2_n500), stack = FALSE)
image_write(g2_combined, "old/g2_violin_combined.png")

# g3组 - J参数变化
g3_combined <- image_append(c(g3_J20, g3_J50, g3_J100), stack = FALSE)
image_write(g3_combined, "old/g3_violin_combined.png")

cat("拼接完成！\n")
cat("- 所有图拼接: old/all_violin_plots_combined.png\n")
cat("- g1组拼接: old/g1_violin_combined.png\n")
cat("- g2组拼接: old/g2_violin_combined.png\n")
cat("- g3组拼接: old/g3_violin_combined.png\n")

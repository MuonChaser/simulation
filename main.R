J=10
n.avg=20
source("gen_data.R")
source("simu.R")
g <- ggplot(bias_data, aes(x = Method, y = Bias, fill = Method)) +
  ylim(-5, 5) +
  geom_boxplot() +
  facet_wrap(~ A) +
  theme_minimal() +
  labs(title = "Bias Distribution by Method and A Level",
       x = "Method",
       y = "Bias") +
  theme(legend.position = "none")

ggsave("j_10_n_20.png", plot = g, width = 8, height = 6)

J=20
n.avg=40
source("gen_data.R")
source("simu.R")
g <- ggplot(bias_data, aes(x = Method, y = Bias, fill = Method)) +
  ylim(-5, 5) +
  geom_boxplot() +
  facet_wrap(~ A) +
  theme_minimal() +
  # 限制 y 轴范围
  labs(title = "Bias Distribution by Method and A Level",
       x = "Method",
       y = "Bias") +
  theme(legend.position = "none")

ggsave("j_20_n_40.png", plot = g, width = 8, height = 6)

J=40
n.avg=60
source("gen_data.R")
source("simu.R")
g <- ggplot(bias_data, aes(x = Method, y = Bias, fill = Method)) +
  ylim(-5, 5) +
  geom_boxplot() +
  facet_wrap(~ A) +
  theme_minimal() +
  labs(title = "Bias Distribution by Method and A Level",
       x = "Method",
       y = "Bias") +
  theme(legend.position = "none")

ggsave("j_40_n_60.png", plot = g, width = 8, height = 6)

data = assign_strata(data)

prop_table <- aggregate(D ~ B + id, data = data, FUN = mean)
names(prop_table) <- c("K", "id", "接受处理比例")
print(prop_table)


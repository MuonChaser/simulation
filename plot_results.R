## plot_results.R
## Creates two publication-quality figures (Liu & Yang style).
##
## Figure 1 (simulation_bias.pdf): 2-row × 3-col facet grid
##   Rows    = mechanism (a=0, a=1)
##   Columns = 3 scenarios
##   x-axis per column: the parameter that varies in that scenario
##     Scenario 1 col: K = 25 / 50 / 100   (n_jk=10, J=10 fixed)
##     Scenario 2 col: n = 100 / 200 / 500  (K=2,  J=10 fixed)
##     Scenario 3 col: J = 20 / 50 / 100   (K=4,  n_jk=20 fixed)
##   Because x is "free per column" in facet_grid(Mechanism ~ scenario),
##   each scenario column shows only its own 3 parameter values.
##
## Figure 2 (simulation_sd.pdf): SD and RMSE vs. varying parameter
##   Same column layout; x-axis shows the scenario-specific parameter.

library(ggplot2)
library(dplyr)
library(tidyr)

## ---- scenario metadata ----
scenarios_meta <- list(
  list(label="g1_K25",  K=25,  n_jk=10,  J=10,  scenario=1, N=2500,
       param_label="K = 25",  param_val=25),
  list(label="g1_K50",  K=50,  n_jk=10,  J=10,  scenario=1, N=5000,
       param_label="K = 50",  param_val=50),
  list(label="g1_K100", K=100, n_jk=10,  J=10,  scenario=1, N=10000,
       param_label="K = 100", param_val=100),
  list(label="g2_n100", K=2,   n_jk=100, J=10,  scenario=2, N=2000,
       param_label="n = 100", param_val=100),
  list(label="g2_n200", K=2,   n_jk=200, J=10,  scenario=2, N=4000,
       param_label="n = 200", param_val=200),
  list(label="g2_n500", K=2,   n_jk=500, J=10,  scenario=2, N=10000,
       param_label="n = 500", param_val=500),
  list(label="g3_J20",  K=4,   n_jk=20,  J=20,  scenario=3, N=1600,
       param_label="J = 20",  param_val=20),
  list(label="g3_J50",  K=4,   n_jk=20,  J=50,  scenario=3, N=4000,
       param_label="J = 50",  param_val=50),
  list(label="g3_J100", K=4,   n_jk=20,  J=100, scenario=3, N=8000,
       param_label="J = 100", param_val=100)
)

## ---- load all RDS files ----
all_data <- lapply(scenarios_meta, function(s) {
  f <- file.path("results", paste0(s$label, "_results.rds"))
  if (!file.exists(f)) return(NULL)
  d <- readRDS(f)
  data.frame(
    bias_strata_a0 = d[, "bias_strata_a0"],
    bias_group_a0  = d[, "bias_group_a0"],
    bias_strata_a1 = d[, "bias_strata_a1"],
    bias_group_a1  = d[, "bias_group_a1"],
    scenario    = paste0("Scenario ", s$scenario),
    param_label = s$param_label,
    param_val   = s$param_val,
    N           = s$N
  )
})
all_df <- do.call(rbind, Filter(Negate(is.null), all_data))

## ---- param_label factor: ordered by param_val WITHIN each scenario ----
## This ensures K=25<50<100, n=100<200<500, J=20<50<100 display left-to-right.
param_order <- all_df %>%
  distinct(scenario, param_label, param_val) %>%
  arrange(scenario, param_val) %>%
  pull(param_label)
## Levels: K=25, K=50, K=100, n=100, n=200, n=500, J=20, J=50, J=100
## (but only the 3 present in each column are shown, others dropped)

## ---- reshape to long format ----
long_df <- all_df %>%
  pivot_longer(
    cols      = starts_with("bias_"),
    names_to  = "key",
    values_to = "bias"
  ) %>%
  mutate(
    Method    = ifelse(grepl("strata", key), "Stratified", "Unstratified"),
    Mechanism = ifelse(grepl("a0$",    key), "italic(a)==0", "italic(a)==1")
  )

long_df$scenario    <- factor(long_df$scenario, levels = paste0("Scenario ", 1:3))
long_df$Mechanism   <- factor(long_df$Mechanism,
                               levels = c("italic(a)==0", "italic(a)==1"))
long_df$Method      <- factor(long_df$Method, levels = c("Stratified", "Unstratified"))
long_df$param_label <- factor(long_df$param_label, levels = param_order)

## ----------------------------------------------------------------
## Figure 1: bias distribution
## Layout: facet_grid(Mechanism ~ scenario)
##   → rows = a∈{0,1}, columns = scenario
##   → x-axis is "free per column" so each scenario column only shows
##     its own 3 parameter values (no staircase)
## ----------------------------------------------------------------
fig1 <- ggplot(long_df, aes(x = param_label, y = bias, fill = Method, colour = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  geom_boxplot(
    outlier.size  = 0.4,
    outlier.alpha = 0.5,
    width         = 0.55,
    position      = position_dodge(0.70),
    linewidth     = 0.35
  ) +
  ## Swap rows/cols: Mechanism as rows, scenario as columns.
  ## scales="free_x": x-axis is independent per column (each scenario).
  ## drop=TRUE in scale_x_discrete: unused factor levels not shown.
  facet_grid(
    Mechanism ~ scenario,
    scales   = "free_x",
    space    = "free_x",
    labeller = labeller(
      Mechanism = label_parsed,
      scenario  = label_value
    )
  ) +
  scale_x_discrete(drop = TRUE) +
  scale_fill_manual(values = c("Stratified" = "white", "Unstratified" = "grey70")) +
  scale_colour_manual(values = c("Stratified" = "black", "Unstratified" = "black")) +
  labs(x = NULL, y = "Bias", fill = NULL, colour = NULL) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text.x       = element_text(size = 9),
    strip.text.y       = element_text(size = 9),
    axis.text.x        = element_text(size = 8.5),
    legend.position    = "none",
    plot.margin        = margin(4, 6, 2, 4)
  )

ggsave("paper/simulation_bias.pdf", fig1, width = 8.5, height = 5.5, device = "pdf")
cat("Saved: paper/simulation_bias.pdf\n")


## ----------------------------------------------------------------
## Figure 3: violin plot (same layout as Figure 1)
## ----------------------------------------------------------------
fig3 <- ggplot(long_df, aes(x = param_label, y = bias, fill = Method, colour = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  geom_violin(
    position  = position_dodge(0.70),
    width     = 0.65,
    linewidth = 0.35,
    trim      = TRUE
  ) +
  facet_grid(
    Mechanism ~ scenario,
    scales   = "free_x",
    space    = "free_x",
    labeller = labeller(
      Mechanism = label_parsed,
      scenario  = label_value
    )
  ) +
  scale_x_discrete(drop = TRUE) +
  scale_fill_manual(values = c("Stratified" = "white", "Unstratified" = "grey70")) +
  scale_colour_manual(values = c("Stratified" = "black", "Unstratified" = "black")) +
  labs(x = NULL, y = "Bias", fill = NULL, colour = NULL) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text.x       = element_text(size = 9),
    strip.text.y       = element_text(size = 9),
    axis.text.x        = element_text(size = 8.5),
    legend.position    = "none",
    plot.margin        = margin(4, 6, 2, 4)
  )

ggsave("paper/simulation_violin.pdf", fig3, width = 8.5, height = 5.5, device = "pdf")
cat("Saved: paper/simulation_violin.pdf\n")


## ----------------------------------------------------------------
## Figure 2: SD & RMSE vs. varying parameter (line plot)
## Same column layout so x-axis shows scenario-specific parameter.
## ----------------------------------------------------------------
stats_data <- lapply(scenarios_meta, function(s) {
  f <- file.path("results", paste0(s$label, "_results.rds"))
  if (!file.exists(f)) return(NULL)
  d <- readRDS(f)
  data.frame(
    scenario    = paste0("Scenario ", s$scenario),
    param_label = s$param_label,
    param_val   = s$param_val,
    sd_strata_a0   = sd(d[, "bias_strata_a0"]),
    sd_group_a0    = sd(d[, "bias_group_a0"]),
    rmse_strata_a0 = sqrt(mean(d[, "bias_strata_a0"]^2)),
    rmse_group_a0  = sqrt(mean(d[, "bias_group_a0"]^2)),
    sd_strata_a1   = sd(d[, "bias_strata_a1"]),
    sd_group_a1    = sd(d[, "bias_group_a1"]),
    rmse_strata_a1 = sqrt(mean(d[, "bias_strata_a1"]^2)),
    rmse_group_a1  = sqrt(mean(d[, "bias_group_a1"]^2))
  )
})
stats_df <- do.call(rbind, Filter(Negate(is.null), stats_data))

sd_long <- stats_df %>%
  pivot_longer(
    cols          = starts_with("sd_") | starts_with("rmse_"),
    names_to      = c("metric", "method", "A"),
    names_pattern = "(sd|rmse)_(strata|group)_(a[01])",
    values_to     = "value"
  ) %>%
  mutate(
    Method    = ifelse(method == "strata", "Stratified", "Unstratified"),
    Mechanism = ifelse(A == "a0", "italic(a)==0", "italic(a)==1"),
    Metric    = ifelse(metric == "sd", "SD", "RMSE")
  )

sd_long$scenario    <- factor(sd_long$scenario, levels = paste0("Scenario ", 1:3))
sd_long$param_label <- factor(sd_long$param_label, levels = param_order)
sd_long$Method      <- factor(sd_long$Method, levels = c("Stratified", "Unstratified"))
sd_long$Mechanism   <- factor(sd_long$Mechanism,
                               levels = c("italic(a)==0", "italic(a)==1"))
sd_long$Metric      <- factor(sd_long$Metric, levels = c("SD", "RMSE"))

fig2 <- ggplot(sd_long,
               aes(x = param_label, y = value, group = Method,
                   colour = Method, linetype = Method, shape = Method)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  ## Rows: Metric (SD / RMSE), then within mechanism group
  ## Columns: scenario (each with its own x-axis parameter)
  facet_grid(
    Metric ~ scenario + Mechanism,
    scales   = "free",
    space    = "free_x",
    labeller = labeller(
      Mechanism = label_parsed,
      scenario  = label_value,
      Metric    = label_value
    )
  ) +
  scale_x_discrete(drop = TRUE) +
  scale_colour_manual(values = c("Stratified" = "black", "Unstratified" = "grey50")) +
  scale_linetype_manual(values = c("Stratified" = "solid", "Unstratified" = "dashed")) +
  scale_shape_manual(values = c("Stratified" = 16, "Unstratified" = 17)) +
  labs(x = NULL, y = NULL, colour = NULL, linetype = NULL, shape = NULL) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text.y       = element_text(size = 9, angle = 0, hjust = 0),
    strip.text.x       = element_text(size = 8),
    axis.text.x        = element_text(size = 7.5),
    legend.position    = "none",
    plot.margin        = margin(4, 6, 2, 4)
  )

ggsave("paper/simulation_sd.pdf", fig2, width = 10, height = 4.5, device = "pdf")
cat("Saved: paper/simulation_sd.pdf\n")

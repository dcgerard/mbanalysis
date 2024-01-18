library(tidyverse)
library(scales)

gdf <- read_csv("./output/sims/gsims.csv")
gdf$rd <- Inf
gldf <- read_csv("./output/sims/glsims.csv")

df <- bind_rows(gdf, gldf)
df |>
  mutate(xi1_c = recode(xi1_c,
                        "low" = "a",
                        "med" = "b",
                        "high" = "c"),
         xi2_c = recode(xi2_c,
                        "low" = "a",
                        "med" = "b",
                        "high" = "c")) |>
  mutate(pgeno = paste0("(", ell1, ",", ell2, ")"),
         Condition = paste0("list(alpha==", alpha_c, ",xi[1]==", xi1_c, ",xi[2]==",xi2_c, ")"),
         Condition = str_remove(Condition, ",xi\\[1\\]==NA"),
         Condition = str_remove(Condition, ",xi\\[2\\]==NA")) ->
  df

## new LRT QQ-plots ----
df |>
  filter(is.infinite(rd)) |>
  ggplot(aes(sample = p_lrt, color = Condition)) +
  facet_grid(pgeno ~ n) +
  geom_qq(distribution = stats::qunif) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        legend.text.align = 0) +
  scale_color_discrete(labels = parse_format()) +
  guides(color = guide_legend(ncol = 1)) +
  xlab("Uniform Quantiles") +
  ylab("Sample P-values") ->
  pl

ggsave(
  filename = "./output/sims/plots/qq_lrt_g.pdf",
  plot = pl,
  height = 6,
  width = 6,
  family = "Times")

df |>
  filter(!is.infinite(rd)) |>
  ggplot(aes(sample = p_lrt, color = Condition)) +
  facet_grid(pgeno ~ n) +
  geom_qq(distribution = stats::qunif) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        legend.text.align = 0) +
  scale_color_discrete(labels = parse_format()) +
  guides(color = guide_legend(ncol = 1)) +
  xlab("Uniform Quantiles") +
  ylab("Sample P-values") ->
  pl

ggsave(
  filename = "./output/sims/plots/qq_lrt_gl.pdf",
  plot = pl,
  height = 6,
  width = 6,
  family = "Times")

## Chi-squared QQ-plots ----
df |>
  filter(is.infinite(rd)) |>
  ggplot(aes(sample = p_chisq, color = Condition)) +
  facet_grid(pgeno ~ n) +
  geom_qq(distribution = stats::qunif) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        legend.text.align = 0) +
  scale_color_discrete(labels = parse_format()) +
  guides(color = guide_legend(ncol = 1)) +
  xlab("Uniform Quantiles") +
  ylab("Sample P-values") ->
  pl

ggsave(
  filename = "./output/sims/plots/qq_chisq_g.pdf",
  plot = pl,
  height = 6,
  width = 6,
  family = "Times")

df |>
  filter(!is.infinite(rd)) |>
  ggplot(aes(sample = p_chisq, color = Condition)) +
  facet_grid(pgeno ~ n) +
  geom_qq(distribution = stats::qunif) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        legend.text.align = 0) +
  scale_color_discrete(labels = parse_format()) +
  guides(color = guide_legend(ncol = 1)) +
  xlab("Uniform Quantiles") +
  ylab("Sample P-values") ->
  pl

ggsave(
  filename = "./output/sims/plots/qq_chisq_gl.pdf",
  plot = pl,
  height = 6,
  width = 6,
  family = "Times")

## polymapR QQ-plots ----
df |>
  filter(is.infinite(rd)) |>
  ggplot(aes(sample = p_polymapr, color = Condition)) +
  facet_grid(pgeno ~ n) +
  geom_qq(distribution = stats::qunif) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        legend.text.align = 0) +
  scale_color_discrete(labels = parse_format()) +
  guides(color = guide_legend(ncol = 1)) +
  xlab("Uniform Quantiles") +
  ylab("Sample P-values") ->
  pl

ggsave(
  filename = "./output/sims/plots/qq_polymapr_g.pdf",
  plot = pl,
  height = 6,
  width = 6,
  family = "Times")

df |>
  filter(!is.infinite(rd)) |>
  ggplot(aes(sample = p_polymapr, color = Condition)) +
  facet_grid(pgeno ~ n) +
  geom_qq(distribution = stats::qunif) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        legend.text.align = 0) +
  scale_color_discrete(labels = parse_format()) +
  guides(color = guide_legend(ncol = 1)) +
  xlab("Uniform Quantiles") +
  ylab("Sample P-values") ->
  pl

ggsave(
  filename = "./output/sims/plots/qq_polymapr_gl.pdf",
  plot = pl,
  height = 6,
  width = 6,
  family = "Times")

## Bayesian boxplots ----
df |>
  filter(is.infinite(rd)) |>
  mutate(`Sample Size` = factor(n)) |>
  ggplot(aes(y = lbf, x = Condition, color = `Sample Size`)) +
  facet_grid(pgeno ~ .) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(labels = parse_format()) +
  ylab("Log Bayes Factor") ->
  pl

ggsave(
  filename = "./output/sims/plots/box_lbf_g.pdf",
  plot = pl,
  height = 6,
  width = 6,
  family = "Times")

df |>
  filter(!is.infinite(rd)) |>
  mutate(`Sample Size` = factor(n)) |>
  ggplot(aes(y = lbf, x = Condition, color = `Sample Size`)) +
  facet_grid(pgeno ~ .) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(labels = parse_format()) +
  ylab("Log Bayes Factor") ->
  pl

ggsave(
  filename = "./output/sims/plots/box_lbf_gl.pdf",
  plot = pl,
  height = 6,
  width = 6,
  family = "Times")

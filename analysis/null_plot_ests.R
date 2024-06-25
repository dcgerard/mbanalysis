library(tidyverse)
library(scales)
pal <- palette.colors(n = 2, palette = "Okabe-Ito", recycle = FALSE, names = FALSE)

gdf <- read_csv("./output/sims/gsims.csv")
gdf$rd <- Inf
gldf <- read_csv("./output/sims/glsims.csv")

df <- bind_rows(gdf, gldf)
df |>
  filter(ell1 != 2 & ell2 != 2) |>
  mutate(pgeno = paste0("(", ell1, ",", ell2, ")")) ->
  df

df |>
  mutate(n = factor(n)) |>
  select(n, alpha_true = alpha, alpha_c, pgeno, alpha_lrt, alpha_bayes) |>
  pivot_longer(cols = c("alpha_lrt", "alpha_bayes"), names_to = "Method", values_to = "alpha") |>
  mutate(Method = recode(Method,
                         "alpha_bayes" = "Bayes",
                         "alpha_lrt" = "MLE")) ->
  dfsub

dfsub |>
  select(alpha_true, alpha_c, pgeno) |>
  distinct() ->
  df_alpha

dfsub |>
  ggplot(aes(x = n, y = alpha, color = Method)) +
  geom_boxplot() +
  facet_grid(alpha_c ~ pgeno) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_manual(values = pal) +
  xlab("Sample Size") +
  ylab(expression(hat(alpha))) +
  geom_hline(data = df_alpha, mapping = aes(yintercept = alpha_true), lty = 2) ->
  pl

ggsave(
  filename = "./output/sims/plots/alpha_ests.pdf",
  plot = pl,
  height = 6,
  width = 6,
  family = "Times")

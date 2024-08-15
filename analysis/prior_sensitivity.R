library(tidyverse)
library(scales)
pal <- palette.colors(n = 7, palette = "Okabe-Ito", recycle = FALSE, names = FALSE)[-c(5, 6)]

## Null Plots
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

df |>
  select(n, ell1, ell2, alpha_c, xi1_c, xi2_c, Condition, rd, pgeno, starts_with("lbf")) |>
  pivot_longer(cols = starts_with("lbf"), names_to = "Prior", values_to = "lbf") |>
  mutate(Prior = recode(Prior,
                        "lbf" = "Default",
                        "lbf_1" = "(1/2,1/2),(1,2)",
                        "lbf_2" = "(1/2,1/2),(1/3,2/3)",
                        "lbf_3" = "(2,2),(1,2)",
                        "lbf_4" = "(2,2),(1/3,2/3)")) ->
  bdf

bdf |>
  filter(is.infinite(rd), n == 20) |>
  ggplot(aes(y = lbf, x = Condition, color = Prior)) +
  facet_grid(pgeno ~ .) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(labels = parse_format()) +
  ylab("Log Bayes Factor") +
  scale_color_manual(values = pal) ->
  pl

ggsave(
  filename = "./output/sims/plots/box_lbf_p_g_20.pdf",
  plot = pl,
  height = 6,
  width = 6,
  family = "Times")

bdf |>
  filter(is.infinite(rd), n == 200) |>
  ggplot(aes(y = lbf, x = Condition, color = Prior)) +
  facet_grid(pgeno ~ .) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(labels = parse_format()) +
  ylab("Log Bayes Factor") +
  scale_color_manual(values = pal) ->
  pl

ggsave(
  filename = "./output/sims/plots/box_lbf_p_g_200.pdf",
  plot = pl,
  height = 6,
  width = 6,
  family = "Times")

bdf |>
  filter(!is.infinite(rd), n == 20) |>
  ggplot(aes(y = lbf, x = Condition, color = Prior)) +
  facet_grid(pgeno ~ .) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(labels = parse_format()) +
  ylab("Log Bayes Factor") +
  scale_color_manual(values = pal) ->
  pl

ggsave(
  filename = "./output/sims/plots/box_lbf_p_gl_20.pdf",
  plot = pl,
  height = 6,
  width = 6,
  family = "Times")

bdf |>
  filter(!is.infinite(rd), n == 200) |>
  ggplot(aes(y = lbf, x = Condition, color = Prior)) +
  facet_grid(pgeno ~ .) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(labels = parse_format()) +
  ylab("Log Bayes Factor") +
  scale_color_manual(values = pal) ->
  pl

ggsave(
  filename = "./output/sims/plots/box_lbf_p_gl_200.pdf",
  plot = pl,
  height = 6,
  width = 6,
  family = "Times")

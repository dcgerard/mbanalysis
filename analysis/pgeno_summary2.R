set.seed(1)
library(updog)
library(tidyverse)
library(plotROC)

## Alt sims
gdf_alt <- read_csv("./output/sims/g_altsims.csv")
gldf_alt <- read_csv("./output/sims/gl_altsims.csv")
gdf_alt$rd <- Inf
gdf_alt$class <- 0
gldf_alt$rd <- 10
gldf_alt$class <- 0

## Null Plots
gdf <- read_csv("./output/sims/gsims.csv")
gdf$rd <- Inf
gdf$class <- 1
gldf <- read_csv("./output/sims/glsims.csv")
gldf$class <- 1

## Subsample Alt Sims to look like table proportions
altprop <- 0.04
gldf_alt |>
  sample_n(size = round(altprop * nrow(gldf))) ->
  gldf_alt_sub

gdf_alt |>
  sample_n(size = round(altprop * nrow(gdf))) ->
  gdf_alt_sub

# Combine
df <- bind_rows(gdf, gldf, gdf_alt_sub, gldf_alt_sub)

pal <- palette.colors(n = 4, palette = "Okabe-Ito")
df |>
  select(n, rd, class, LRT = p_lrt, Chisq = p_chisq, polymapR = p_polymapr, Bayes = lbf) |>
  pivot_longer(cols = c("LRT", "Chisq", "polymapR", "Bayes"), names_to = "Method") |>
  mutate(n = paste0("n=", n),
         rd = paste0("rd=", rd)) |>
  ggplot(aes(m = value, d = class, color = Method)) +
  geom_roc(labels = FALSE) +
  facet_grid(n ~ rd) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  scale_color_manual(values = pal) +
  xlim(0, 0.1) ->
  pl

ggsave(
  filename = "./output/roc/roc_plot_2.pdf",
  plot = pl,
  height = 4,
  width = 6,
  family = "Times")

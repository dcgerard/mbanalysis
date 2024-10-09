library(updog)
library(tidyverse)
library(plotROC)

## Get Empirical Proportions
bout <- readRDS(file = "./output/blue/bluefits.RDS")
geno <- format_multidog(bout, "geno")
gsub <- geno[, colnames(geno) %in% c("indigocrisp", "sweetcrisp")]

gsub[gsub == 3] <- 1
gsub[gsub == 4] <- 0

tab <- table(gsub[, "indigocrisp"], gsub[, "sweetcrisp"])
dimnames(tab) <- list("Parent1" = c("nullplex", "simplex", "duplex"),
                      "Parent2" = c("nullplex", "simplex", "duplex"))
data.frame(ell1 = c(0, 0, 1, 1, 2),
           ell2 = c(1, 2, 1, 2, 2),
           num = c(tab[1, 2] + tab[2, 1],
                   tab[1, 3] + tab[3, 1],
                   tab[2, 2],
                   tab[2, 3] + tab[3, 2],
                   tab[3, 3])) |>
  mutate(prop = num / sum(num)) ->
  tabdf

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

## Subsample Null Sims to look like table proportions
gdf |>
  group_by(ell1, ell2) |>
  summarize(n = n()) |>
  ungroup() |>
  mutate(prop = n / sum(n)) ->
simdf

tabdf |>
  mutate(nsamp = round(prop * 4108)) ->
  tabdf

gdf |>
  group_by(ell1, ell2) |>
  nest() |>
  left_join(tabdf, by = join_by(ell1, ell2)) |>
  mutate(subdata = map2(.x = data, .y = nsamp, .f = \(x, y) slice_sample(.data = x, n = y))) |>
  select(ell1, ell2, subdata) |>
  unnest(subdata) ->
  gdf_sub

gldf |>
  group_by(ell1, ell2) |>
  nest() |>
  left_join(tabdf, by = join_by(ell1, ell2)) |>
  mutate(subdata = map2(.x = data, .y = nsamp, .f = \(x, y) slice_sample(.data = x, n = y))) |>
  select(ell1, ell2, subdata) |>
  unnest(subdata) ->
  gldf_sub

# Combine
df <- bind_rows(gdf_sub, gldf_sub, gdf_alt, gldf_alt)

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
  filename = "./output/roc/roc_plot_1.pdf",
  plot = pl,
  height = 4,
  width = 6,
  family = "Times")

## Generate a simplex grid and plot with ggtern
library(tidyverse)
library(ggtern)
library(latex2exp)

## Alternative data frame
np <- 6
df_alt <- expand.grid(
  p0 = seq(0, 1, length.out = np),
  p1 = seq(0, 1, length.out = np),
  p2 = seq(0, 1, length.out = np)) |>
  filter(p0 + p1 + p2 == 1) |>
  filter(!(p0 == p2 & p1 > 0.5),
         !(p0 == 0 & p1 == 0 & p2 == 1),
         !(p0 == 0 & p1 == 1 & p2 == 0),
         !(p0 == 1 & p1 == 0 & p2 == 0))

## Null Data Frames
df0 <- matrix(menbayes::pvec_tet_2(alpha = 0, xi = 1/3, ell = 0), ncol = 3) |>
  as.data.frame()
names(df0) <- c("p0", "p1", "p2")

aseq <- seq(0, 1/6, length.out = 100)
df1 <- matrix(menbayes::pvec_tet_2(alpha = aseq, xi = 1/3, ell = 1), ncol = 3) |>
  as.data.frame()
names(df1) <- c("p0", "p1", "p2")

xiseq <- seq(0, 1, length.out = 100)
df2 <- matrix(menbayes::pvec_tet_2(alpha = 0, xi = xiseq, ell = 2), ncol = 3) |>
  as.data.frame()
names(df2) <- c("p0", "p1", "p2")

aseq <- seq(0, 1/6, length.out = 100)
df3 <- matrix(menbayes::pvec_tet_2(alpha = aseq, xi = 1/3, ell = 3), ncol = 3) |>
  as.data.frame()
names(df3) <- c("p0", "p1", "p2")

df4 <- matrix(menbayes::pvec_tet_2(alpha = 0, xi = 1/3, ell = 4), ncol = 3) |>
  as.data.frame()
names(df4) <- c("p0", "p1", "p2")

# Make the plot
ggtern(data = df_alt, mapping = aes(x = p0, y = p1, z = p2)) +
  geom_point(size = 3) +
  geom_point(data = df0, color = "blue", size = 7) +
  geom_line(data = df1, color = "blue", lwd = 2) +
  geom_line(data = df2, color = "blue", lwd = 2) +
  geom_line(data = df3, color = "blue", lwd = 2) +
  geom_point(data = df4, color = "blue", size = 7) +
  theme_bw() +
  scale_L_continuous(name = TeX("$p_0$"), labels = labels_tern(factor = 1)) +
  scale_T_continuous(name = TeX("$p_1$"), labels = labels_tern(factor = 1)) +
  scale_R_continuous(name = TeX("$p_2$"), labels = labels_tern(factor = 1)) ->
  pl

ggsave(filename = "./output/hyp/alt_studied.pdf", plot = pl, height = 4, width = 4)

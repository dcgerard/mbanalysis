## Plot various hypotheses in 2-simplex
library(tidyverse)
library(ggtern)
library(menbayes)
library(latex2exp)

## ell = 2
df2_1 <- matrix(menbayes::pvec_tet_2(alpha = 0, xi = 1/3, ell = 2), ncol = 3) |>
  as.data.frame()
names(df2_1) <- c("p0", "p1", "p2")

aseq <- seq(0, 1/6, length.out = 100)
df2_2 <- matrix(menbayes::pvec_tet_2(alpha = aseq, xi = 1/3, ell = 2), ncol = 3) |>
  as.data.frame()
names(df2_2) <- c("p0", "p1", "p2")

xiseq <- seq(0, 1, length.out = 100)
df2_3 <- matrix(menbayes::pvec_tet_2(alpha = 0, xi = xiseq, ell = 2), ncol = 3) |>
  as.data.frame()
names(df2_3) <- c("p0", "p1", "p2")


pl <- ggtern(mapping = aes(x = p0, y = p1, z = p2)) +
  geom_line(data = df2_3, color = "black", lwd = 1.5) +
  geom_line(data = df2_2, color = "red", lwd = 1.5) +
  geom_point(data = df2_1, color = "blue", cex = 3) +
  theme_bw() +
  scale_L_continuous(name = TeX("$p_0$"), labels = labels_tern(factor = 1)) +
  scale_T_continuous(name = TeX("$p_1$"), labels = labels_tern(factor = 1)) +
  scale_R_continuous(name = TeX("$p_2$"), labels = labels_tern(factor = 1))

# ell = 1
df1_1 <- matrix(menbayes::pvec_tet_2(alpha = 0, xi = 1/3, ell = 1), ncol = 3) |>
  as.data.frame()
names(df1_1) <- c("p0", "p1", "p2")

aseq <- seq(0, 1/6, length.out = 100)
df1_2 <- matrix(menbayes::pvec_tet_2(alpha = aseq, xi = 1/3, ell = 1), ncol = 3) |>
  as.data.frame()
names(df1_2) <- c("p0", "p1", "p2")

pl <- pl +
  geom_line(data = df1_2, color = "red", lwd = 1.5) +
  geom_point(data = df1_1, color = "blue", cex = 3)

# ell = 3
df3_1 <- matrix(menbayes::pvec_tet_2(alpha = 0, xi = 1/3, ell = 3), ncol = 3) |>
  as.data.frame()
names(df3_1) <- c("p0", "p1", "p2")

aseq <- seq(0, 1/6, length.out = 100)
df3_2 <- matrix(menbayes::pvec_tet_2(alpha = aseq, xi = 1/3, ell = 3), ncol = 3) |>
  as.data.frame()
names(df3_2) <- c("p0", "p1", "p2")

pl <- pl +
  geom_line(data = df3_2, color = "red", lwd = 1.5) +
  geom_point(data = df3_1, color = "blue", cex = 3)

# ell = 4
df4_1 <- matrix(menbayes::pvec_tet_2(alpha = 0, xi = 1/3, ell = 4), ncol = 3) |>
  as.data.frame()
names(df4_1) <- c("p0", "p1", "p2")

pl <- pl +
  geom_point(data = df4_1, color = "blue", cex = 5)

# ell = 0
df0_1 <- matrix(menbayes::pvec_tet_2(alpha = 0, xi = 1/3, ell = 0), ncol = 3) |>
  as.data.frame()
names(df0_1) <- c("p0", "p1", "p2")

pl <- pl +
  geom_point(data = df0_1, color = "blue", cex = 5)

ggsave(filename = "./output/hyp/ternary.pdf", plot = pl, height = 4, width = 4, family = "Times")

## Make material for legend
dftrash <- data.frame(x = 1:3, y = 1:3, z = c("No PP, No DR", "DR, No PP", "DR, PP"))
ggplot(dftrash, aes(x = x, y = y, color = z)) +
  geom_line(lwd = 1.5) +
  theme_bw() +
  scale_color_manual(values =  c("red", "black", "blue")) ->
  pl
ggsave(filename = "./output/hyp/line_legend.pdf", plot = pl, height = 4, width = 4)

ggplot(dftrash, aes(x = x, y = y, color = z)) +
  geom_point(cex = 3) +
  theme_bw() +
  scale_color_manual(values =  c("red", "black", "blue")) ->
  pl
ggsave(filename = "./output/hyp/point_legend.pdf", plot = pl, height = 4, width = 4)


tibble(x = 1:5, y = 1:5, text = c("l==0", "l==1", "l==2", "\u2113==3", "\u2113==4")) |>
  ggplot(aes(x = x, y = y, label = text)) +
  geom_text(parse = TRUE) +
  theme_void() ->
  pl
ggsave(filename = "./output/hyp/ell_text.pdf", plot = pl, height = 4, width = 4)

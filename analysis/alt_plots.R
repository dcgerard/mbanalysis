## Alternative plots
library(tidyverse)
library(xtable)

gdf <- read_csv("./output/sims/g_altsims.csv")
gldf <- read_csv("./output/sims/gl_altsims.csv")

gdf$rd <- "Inf"
gldf$rd <- "10"

df <- bind_rows(gdf, gldf)

# Create data frame of alternatives
np <- 6
qmat <- rbind(
  rep(0.2, 5),
  c(1/4, 1/4, 1/4, 1/4, 0),
  c(4/10, 3/10, 2/10, 1/10, 0),
  c(1/10, 2/10, 3/10, 4/10, 0),
  c(1/3, 1/3, 1/3, 0, 0),
  c(0, 1/3, 1/3, 1/3, 0),
  c(3/6, 2/6, 1/6, 0, 0),
  c(0, 3/6, 2/6, 1/6, 0),
  c(1/6, 2/6, 3/6, 0, 0),
  c(0, 1/6, 2/6, 3/6, 0),
  c(3/4, 1/4, 0, 0, 0),
  c(1/4, 3/4, 0, 0, 0),
  c(0, 3/4, 1/4, 0, 0),
  c(0, 1/4, 3/4, 0, 0))
colnames(qmat) <- c("q0", "q1", "q2", "q3", "q4")
as_tibble(qmat) |>
  mutate(name = str_c("alt_", 1:n())) |>
  mutate(q = str_c("(", str_c(round(q0, digits = 2), round(q1, digits = 2), round(q2, digits = 2), round(q3, digits = 2), round(q4, digits = 2), sep = ", "), ")"))

## Analysis of p-values --------------------------------------------------------
## Vary n, rd, alt, create a power vs t1e plot
df |>
  select(n, rd, name = alt, LRT = p_lrt, Chisq = p_chisq, polymapR = p_polymapr) |>
  left_join(df_alt, by = join_by(name)) |>
  mutate(q = case_when(name == "random" ~ "Random",
                       TRUE ~ q)) |>
  select(n, rd, q, LRT, Chisq, polymapR) ->
  df_roc

# p to power function
#' @param p A vector of p-values
#' @param alpha A vector of significance levels
#'
#' @return A vector of powers the same length as alpha
p_to_power <- function(p, alpha) {
  power <- rep(NA_real_, length.out = length(alpha))
  for (i in seq_along(alpha)) {
    power[[i]] <- mean(p < alpha[[i]], na.rm = TRUE)
  }
  return(power)
}

aseq <- seq(1e-5, 0.1, length.out = 100)
pal <- palette.colors(n = 3, palette = "Okabe-Ito")
df_roc |>
  filter(n == 20, rd == "10") |>
  select(-n, -rd) |>
  pivot_longer(cols = c("LRT", "Chisq", "polymapR"), names_to = "Method", values_to = "p_value") |>
  group_by(q, Method) |>
  reframe(alpha = aseq, power = p_to_power(p = p_value, alpha = aseq)) |>
  ggplot(aes(x = alpha, y = power, color = Method)) +
  geom_line() +
  facet_wrap(. ~ q) +
  ylab("Power") +
  xlab("Stated Significance Level") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text = element_text(size = 6)) +
  scale_color_manual(values = pal) ->
  pl

ggsave(
  filename = "./output/sims/alt_plots/alt_power_n20_rd10.pdf",
  plot = pl,
  height = 5,
  width = 6,
  family = "Times")

df_roc |>
  filter(n == 200, rd == "10") |>
  select(-n, -rd) |>
  pivot_longer(cols = c("LRT", "Chisq", "polymapR"), names_to = "Method", values_to = "p_value") |>
  group_by(q, Method) |>
  reframe(alpha = aseq, power = p_to_power(p = p_value, alpha = aseq)) |>
  ggplot(aes(x = alpha, y = power, color = Method)) +
  geom_line() +
  facet_wrap(. ~ q) +
  ylab("Power") +
  xlab("Stated Significance Level") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text = element_text(size = 6)) +
  scale_color_manual(values = pal) ->
  pl

ggsave(
  filename = "./output/sims/alt_plots/alt_power_n200_rd10.pdf",
  plot = pl,
  height = 5,
  width = 6,
  family = "Times")

df_roc |>
  filter(n == 20, rd == "Inf") |>
  select(-n, -rd) |>
  pivot_longer(cols = c("LRT", "Chisq", "polymapR"), names_to = "Method", values_to = "p_value") |>
  group_by(q, Method) |>
  reframe(alpha = aseq, power = p_to_power(p = p_value, alpha = aseq)) |>
  ggplot(aes(x = alpha, y = power, color = Method)) +
  geom_line() +
  facet_wrap(. ~ q) +
  ylab("Power") +
  xlab("Stated Significance Level") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text = element_text(size = 6)) +
  scale_color_manual(values = pal) ->
  pl

ggsave(
  filename = "./output/sims/alt_plots/alt_power_n20_rdInf.pdf",
  plot = pl,
  height = 5,
  width = 6,
  family = "Times")

df_roc |>
  filter(n == 200, rd == "Inf") |>
  select(-n, -rd) |>
  pivot_longer(cols = c("LRT", "Chisq", "polymapR"), names_to = "Method", values_to = "p_value") |>
  group_by(q, Method) |>
  reframe(alpha = aseq, power = p_to_power(p = p_value, alpha = aseq)) |>
  ggplot(aes(x = alpha, y = power, color = Method)) +
  geom_line() +
  facet_wrap(. ~ q) +
  ylab("Power") +
  xlab("Stated Significance Level") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text = element_text(size = 6)) +
  scale_color_manual(values = pal) ->
  pl

ggsave(
  filename = "./output/sims/alt_plots/alt_power_n200_rdInf.pdf",
  plot = pl,
  height = 5,
  width = 6,
  family = "Times")

## Analysis of Bayes factors --------------------------------------------------

pal <- palette.colors(n = 5, palette = "Okabe-Ito")
df |>
  select(n, rd, name = alt, starts_with("lbf")) |>
  mutate(Condition = paste0("n=", n, ",rd=", rd)) |>
  left_join(df_alt, by = join_by(name)) |>
  mutate(q = case_when(name == "unif" ~ "(.2, .2, .2, .2, .2)",
                       name == "random" ~ "Random",
                       TRUE ~ q)) |>
  select(Condition, q, starts_with("lbf")) |>
  pivot_longer(cols = starts_with("lbf"), names_to = "Prior", values_to = "lbf") |>
  mutate(Prior = recode(Prior,
                        "lbf" = "Default",
                        "lbf_1" = "(1/2,1/2),(1,2)",
                        "lbf_2" = "(1/2,1/2),(1/3,2/3)",
                        "lbf_3" = "(2,2),(1,2)",
                        "lbf_4" = "(2,2),(1/3,2/3)")) |>
  ggplot(aes(x = Condition, y = lbf, color = Prior)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(.~q, scales = "free_y") +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text = element_text(size = 6)) +
  scale_color_manual(values = pal) +
  ylab("Log Bayes Factor") ->
  pl

ggsave(
  filename = "./output/sims/alt_plots/alt_lbf_box.pdf",
  plot = pl,
  height = 6,
  width = 8,
  family = "Times")

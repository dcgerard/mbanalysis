library(tidyverse)
bdf <- read_csv("./output/blue/blue_df.csv")
bdf |>
  select(snp, p1, p2, lrt_alpha, lrt_xi1, lrt_xi2) |>
  separate(col = "snp", into = c("chrom", "pos"), sep = "_") |>
  mutate(chrom = parse_number(chrom), pos = parse_number(pos)) ->
  bdf

## Hypothesis Tests, parent 1
nbreaks <- 20
z <- stats::qnorm(1 - 0.05/2)
bdf |>
  filter(chrom <= 22) |>
  filter(p1 %in% c(1, 3) & p2 %in% c(0, 4)) |>
  group_by(chrom) |>
  mutate(pos_cat = cut(x = pos, breaks = nbreaks, labels = 1:nbreaks)) |>
  mutate(pos = as.numeric(pos_cat)) |>
  mutate(quant = case_when(pos %in% c(1, 2) ~ "A",
                           pos %in% c(19, 20) ~ "C",
                           pos %in% c(9, 10, 11, 12) ~ "B",
                           TRUE ~ NA)) |>
  filter(!is.na(quant)) ->
  bdf_sub

bdf_sub |>
  group_by(chrom) |>
  nest() |>
  mutate(aov = map(data, \(x) aov(lrt_alpha ~ quant, data = x))) |>
  mutate(tukey = map(aov, \(x) TukeyHSD(x = x)$quant[, 4])) |>
  unnest(tukey) |>
  mutate(comparison = c("B-A", "C-A", "C-B")) |>
  filter(comparison %in% c("B-A", "C-B")) |>
  arrange(chrom) |>
  mutate(quant = recode(comparison, "B-A" = "A", "C-B" = "C")) |>
  select(chrom, pvalue = tukey, quant) ->
  df_sig

bdf_sub |>
  group_by(chrom, quant) |>
  summarize(m_alpha = mean(lrt_alpha), se_alpha = sd(lrt_alpha) / sqrt(n())) |>
  mutate(lower = m_alpha - z * se_alpha, upper = m_alpha + z * se_alpha) ->
  df_mean

df_mean |>
  left_join(df_sig, by = join_by(chrom, quant)) |>
  mutate(pvalue = format(pvalue, digits = 2),
         pvalue = str_replace(pvalue, "NA", "")) |>
  ggplot(aes(x = quant, y = m_alpha)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = lower, ymax = upper)) +
  facet_wrap(. ~ chrom, scales = "free_x") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  geom_text(aes(x = quant, y = upper, label = pvalue), nudge_y = 0.02) ->
  pl

ggsave("./output/blue/plots/tukey_p1.pdf", plot = pl, height = 6, width = 6, family = "Times")

## Hypothesis Tests, parent 2
nbreaks <- 20
z <- stats::qnorm(1 - 0.05/2)
bdf |>
  filter(chrom <= 22) |>
  filter(p1 %in% c(0, 4) & p2 %in% c(1, 3)) |>
  group_by(chrom) |>
  mutate(pos_cat = cut(x = pos, breaks = nbreaks, labels = 1:nbreaks)) |>
  mutate(pos = as.numeric(pos_cat)) |>
  mutate(quant = case_when(pos %in% c(1, 2) ~ "A",
                           pos %in% c(19, 20) ~ "C",
                           pos %in% c(9, 10, 11, 12) ~ "B",
                           TRUE ~ NA)) |>
  filter(!is.na(quant)) ->
  bdf_sub

bdf_sub |>
  group_by(chrom) |>
  nest() |>
  mutate(aov = map(data, \(x) aov(lrt_alpha ~ quant, data = x))) |>
  mutate(tukey = map(aov, \(x) TukeyHSD(x = x)$quant[, 4])) |>
  unnest(tukey) |>
  mutate(comparison = c("B-A", "C-A", "C-B")) |>
  filter(comparison %in% c("B-A", "C-B")) |>
  arrange(chrom) |>
  mutate(quant = recode(comparison, "B-A" = "A", "C-B" = "C")) |>
  select(chrom, pvalue = tukey, quant) ->
  df_sig

bdf_sub |>
  group_by(chrom, quant) |>
  summarize(m_alpha = mean(lrt_alpha), se_alpha = sd(lrt_alpha) / sqrt(n())) |>
  mutate(lower = m_alpha - z * se_alpha, upper = m_alpha + z * se_alpha) ->
  df_mean

df_mean |>
  left_join(df_sig, by = join_by(chrom, quant)) |>
  mutate(pvalue = format(pvalue, digits = 2),
         pvalue = str_replace(pvalue, "NA", "")) |>
  ggplot(aes(x = quant, y = m_alpha)) +
  geom_point() +
  geom_errorbar(mapping = aes(ymin = lower, ymax = upper)) +
  facet_wrap(. ~ chrom, scales = "free_x") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  geom_text(aes(x = quant, y = upper, label = pvalue), nudge_y = 0.02) ->
  pl

ggsave("./output/blue/plots/tukey_p2.pdf", plot = pl, height = 6, width = 6, family = "Times")

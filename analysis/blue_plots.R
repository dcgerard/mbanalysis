library(tidyverse)

bdf <- read_csv("./output/blue/blue_df.csv")

bdf |>
  mutate(Chromosome = parse_number(str_extract(snp, "^\\d+\\_")),
         Position = parse_number(str_extract(snp, "\\_\\d+$"))) |>
  filter(Chromosome <= 22) ->
  bdf

bdf |>
  mutate(Chromosome = factor(Chromosome)) |>
  ggplot(aes(x = Chromosome, y = lbf)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 2, col = 2) +
  theme_bw() +
  ylab("Log Bayes Factor")

bdf |>
  mutate(Chromosome = factor(Chromosome)) |>
  ggplot(aes(sample = lrt_pvalue, color = Chromosome)) +
  geom_qq(distribution = stats::qunif) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 0, lty = 2, col = 2) +
  theme_bw() +
  ylab("Log Bayes Factor")

bdf |>
  mutate(Chromosome = factor(Chromosome)) |>
  ggplot(aes(sample = polymapr_pvalue, color = Chromosome)) +
  geom_qq(distribution = stats::qunif) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 0, lty = 2, col = 2) +
  theme_bw() +
  ylab("Log Bayes Factor")

## alpha estimates
bdf |>
  filter((p1 %in% c(1, 3) & p2 %in% c(0, 4)) |
           (p1 %in% c(0, 4) & p2 %in% c(1, 3))) |>
  ggplot(aes(x = Position, y = lrt_alpha)) +
  facet_wrap(.~Chromosome) +
  geom_point() +
  geom_smooth()



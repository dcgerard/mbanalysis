library(tidyverse)
library(xtable)

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

# Combine
df <- bind_rows(gdf, gldf, gdf_alt, gldf_alt)

# look at just rd = Inf, n = 200
df |>
  filter(is.infinite(rd), n == 200) |>
  select(ell1, ell2, class, p_lrt, p_polymapr) ->
  df_sub

t1e <- 0.025
df_sub |>
  arrange(p_lrt) |>
  mutate(t1e_lrt = cumsum(class) / sum(class)) |>
  mutate(result_lrt = t1e_lrt <= t1e) |>
  arrange(p_polymapr) |>
  mutate(t1e_polymapr = cumsum(class) / sum(class)) |>
  mutate(result_polymapr = t1e_polymapr <= t1e) |>
  filter(class == 1) |>
  group_by(ell1, ell2) |>
  summarize(nreject_lrt = sum(result_lrt),
            nreject_polymapr = sum(result_polymapr),
            n = n()) ->
  df_sum

df_sum$t1e_lrt <- NA_real_
df_sum$lower_lrt <- NA_real_
df_sum$upper_lrt <- NA_real_

df_sum$t1e_polymapr <- NA_real_
df_sum$lower_polymapr <- NA_real_
df_sum$upper_polymapr <- NA_real_
for (i in seq_len(nrow(df_sum))) {
  bout_lrt <- binom.test(x = df_sum$nreject_lrt[[i]], n = df_sum$n[[i]])
  bout_polymapr <- binom.test(x = df_sum$nreject_polymapr[[i]], n = df_sum$n[[i]])

  df_sum$t1e_lrt[[i]] <- bout_lrt$estimate
  df_sum$lower_lrt[[i]] <- bout_lrt$conf.int[[1]]
  df_sum$upper_lrt[[i]] <- bout_lrt$conf.int[[2]]

  df_sum$t1e_polymapr[[i]] <- bout_polymapr$estimate
  df_sum$lower_polymapr[[i]] <- bout_polymapr$conf.int[[1]]
  df_sum$upper_polymapr[[i]] <- bout_polymapr$conf.int[[2]]
}

df_sum |>
  select(-n, -nreject_lrt, -nreject_polymapr) |>
  xtable(digits = 3) |>
  print(include.rownames = FALSE)

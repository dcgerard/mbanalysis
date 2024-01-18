## Alternative plots
library(tidyverse)
library(xtable)

gdf <- read_csv("./output/sims/g_altsims.csv")
gldf <- read_csv("./output/sims/gl_altsims.csv")

gdf$rd <- "Inf"
gldf$rd <- "10"

df <- bind_rows(gdf, gldf)

df |>
  mutate(Condition = paste0("n=", n, ",rd=", rd)) |>
  select(Condition, LRT = p_lrt, Chisq = p_chisq, polymapR = p_polymapr) |>
  pivot_longer(cols = c("LRT", "Chisq", "polymapR"), names_to = "Method", values_to = "P-value") |>
  ggplot(aes(x = Method, y = `P-value`, color = `Condition`)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() ->
  pl

ggsave(
  filename = "./output/sims/plots/alt_p_box.pdf",
  plot = pl,
  height = 3,
  width = 4,
  family = "Times")

df |>
  mutate(Condition = paste0("n=", n, ",rd=", rd)) |>
  select(Condition, LRT = p_lrt, Chisq = p_chisq, polymapR = p_polymapr) |>
  pivot_longer(cols = c("LRT", "Chisq", "polymapR"), names_to = "Method", values_to = "P-value") |>
  group_by(Condition, Method) |>
  summarize(`.05` = mean(`P-value` < 0.05),
            `.01` = mean(`P-value` < 0.01),
            `.001` = mean(`P-value` < 0.001),
            `.0001` = mean(`P-value` < 0.0001),
            `.00001` = mean(`P-value` < 0.00001)) |>
  xtable(caption = "Power for different methods (Method) at different samples sizes and read-depth (Condition) at different significance levels (.05 through .00001). The methods considered are the standard chi-squared test (Chisq), the likelihood ratio test of Section \\ref{sec:lrt} (LRT), and the \\textsf{polymapR} test of Section \\ref{sec:polymapr}.",
         label = "tab:t1e") |>
  print(file = "./output/sims/plots/t1e.tex",
        include.rownames = FALSE)



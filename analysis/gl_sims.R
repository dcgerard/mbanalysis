library(tidyverse)
library(menbayes)
library(doFuture)
library(doRNG)

## Set up parallelization ----
#  R CMD BATCH '--args nc=8' sims.R
registerDoFuture()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nc <- 1
} else {
  eval(parse(text = args[[1]]))
}
if (nc == 1) {
  plan("sequential")
} else {
  plan("multisession", workers = nc)
}

## Set up parameter settings -----
drbound <- 1/6
pardf <- expand_grid(
  seed = 1:200,
  n = c(25, 100),
  ell1 = 0:2,
  ell2 = 1:2,
  alpha_c = c("0", "1/12", "1/6"),
  xi1_c = c("low", "med", "high"),
  xi2_c = c("low", "med", "high"),
  rd = 10) |>
  filter(!(xi1_c == "high" & xi2_c == "med")) |>
  filter(!(xi1_c == "high" & xi2_c == "low")) |>
  filter(!(xi1_c == "med" & xi2_c == "low")) |>
  filter(!(ell1 == 2 & ell2 == 1)) |>
  filter(!(ell1 == 0 & xi1_c == "high")) |>
  filter(!(ell1 == 1 & xi1_c == "high")) |>
  filter(!(ell2 == 1 & xi2_c == "low")) |>
  filter(!(ell2 == 1 & xi2_c == "high")) |>
  filter(!(alpha_c == "1/6" & xi1_c == "low")) |>
  filter(!(alpha_c == "1/6" & xi1_c == "high")) |>
  filter(!(alpha_c == "1/6" & xi2_c == "low")) |>
  filter(!(alpha_c == "1/6" & xi2_c == "high")) |>
  mutate(alpha = case_when(alpha_c == "0" ~ 0,
                           alpha_c == "1/12" ~ 1/12,
                           alpha_c == "1/6" ~ 1/6)) |>
  mutate(xil = (1/3) * (alpha / (1 - alpha)) * ((1 - drbound) / drbound),
         xiu = 1 - 2 * xil,
         xi1 = case_when(xi1_c == "low" ~ xil,
                         xi1_c == "med" ~ 1/3,
                         xi1_c == "high" ~ xiu),
         xi2 = case_when(xi2_c == "low" ~ xil,
                         xi2_c == "med" ~ 1/3,
                         xi2_c == "high" ~ xiu),
         xi1 = if_else(ell1 != 2, NA_real_, xi1),
         xi2 = if_else(ell2 != 2, NA_real_, xi2),
         xi1_c = if_else(ell1 != 2, NA_character_, xi1_c),
         xi2_c = if_else(ell2 != 2, NA_character_, xi2_c)
  )

## new lrt params
pardf$p_lrt <- NA_real_
pardf$stat_lrt <- NA_real_
pardf$df_lrt <- NA_real_
pardf$alpha_lrt <- NA_real_
pardf$xi1_lrt <- NA_real_
pardf$xi2_lrt <- NA_real_

## lrt_ndr_npp params
pardf$p_lrtnn <- NA_real_
pardf$stat_lrtnn <- NA_real_
pardf$df_lrtnn <- NA_real_

## old chisq params
pardf$p_chisq <- NA_real_
pardf$df_chisq <- NA_real_
pardf$stat_chisq <- NA_real_

## old polymapr params
pardf$p_polymapr <- NA_real_

## new bayes params
pardf$lbf <- NA_real_
pardf$alpha_bayes <- NA_real_
pardf$xi1_bayes <- NA_real_
pardf$xi2_bayes <- NA_real_

## Run simulations
outdf <- foreach(
  i = seq_len(nrow(pardf)),
  .combine = rbind,
  .export = c("pardf")) %dorng% {

    gl <- simf1gl(
      n = pardf$n[[i]],
      g1 = pardf$ell1[[i]],
      g2 = pardf$ell2[[i]],
      rd = pardf$rd[[i]],
      alpha = pardf$alpha[[i]],
      xi1 = pardf$xi1[[i]],
      xi2 = pardf$xi2[[i]])

      ## new lrt ----
      tout <- lrt_men_gl4(gl = gl, g1 = pardf$ell1[[i]], g2 = pardf$ell2[[i]])
      pardf$p_lrt[[i]] <- tout$p_value
      pardf$stat_lrt[[i]] <- tout$statistic
      pardf$df_lrt[[i]] <- tout$df
      pardf$alpha_lrt[[i]] <- tout$alpha
      pardf$xi1_lrt[[i]] <- tout$xi1
      pardf$xi2_lrt[[i]] <- tout$xi2

      ## Basic LRT ----
      tout <- lrt_men_gl4(gl = gl, g1 = pardf$ell1[[i]], g2 = pardf$ell2[[i]], pp = FALSE, dr = FALSE, alpha = 0, xi1 = 1/3, xi2 = 1/3)
      pardf$p_lrtnn[[i]] <- tout$p_value
      pardf$stat_lrtnn[[i]] <- tout$statistic
      pardf$df_lrtnn[[i]] <- tout$df

      ## old chisq ----
      suppressWarnings( ## because of small counts
        nout <- chisq_gl4(gl = gl, g1 = pardf$ell1[[i]], g2 = pardf$ell2[[i]])
      )
      pardf$p_chisq[[i]] <- nout$p_value
      pardf$df_chisq[[i]] <- nout$df
      pardf$stat_chisq[[i]] <- nout$statistic

      ## old polymapr ----
      gp <- exp(gl - apply(X = gl, MARGIN = 1, FUN = updog::log_sum_exp))
      pout <- polymapr_test(x = gp, g1 = pardf$ell1[[i]], g2 = pardf$ell2[[i]], type = "menbayes")
      pardf$p_polymapr[[i]] <- pout$p_value

      ## Bayes test ----
      trash <- capture.output(
        bout <- bayes_men_gl4(gl = gl, g1 = pardf$ell1[[i]], g2 = pardf$ell2[[i]], chains = 1)
      )
      pardf$lbf[[i]] <- bout$lbf
      pardf$alpha_bayes[[i]] <- bout$alpha
      pardf$xi1_bayes[[i]] <- bout$xi1
      pardf$xi2_bayes[[i]] <- bout$xi2

      pardf[i, ]
  }

## Unregister workers ----
if (nc > 1) {
  plan(sequential)
}

write_csv(x = outdf, file = "./output/sims/glsims.csv")

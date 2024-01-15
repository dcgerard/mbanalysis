library(tidyverse)
library(menbayes)

pardf <- expand_grid(
  seed = 1:200,
  n = c(25, 100),
  rd = 10)

qvec <- rep(1/5, 5)

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

for (i in seq_len(nrow(pardf))) {
  cat(i, " of ", nrow(pardf), "\n")
  g1 <- 2
  g2 <- 2
  x <- c(stats::rmultinom(n = 1, size = pardf$n[[i]], prob = qvec))
  gl <- hwep::simgl(nvec = x, rdepth = pardf$rd[[i]], ret = "gl")

  ## new lrt ----
  tout <- lrt_men_gl4(gl = gl, g1 = g1, g2 = g2)
  pardf$p_lrt[[i]] <- tout$p_value
  pardf$stat_lrt[[i]] <- tout$statistic
  pardf$df_lrt[[i]] <- tout$df
  pardf$alpha_lrt[[i]] <- tout$alpha
  pardf$xi1_lrt[[i]] <- tout$xi1
  pardf$xi2_lrt[[i]] <- tout$xi2

  ## Basic LRT ----
  tout <- lrt_men_gl4(gl = gl, g1 = g1, g2 = g2, pp = FALSE, dr = FALSE, alpha = 0, xi1 = 1/3, xi2 = 1/3)
  pardf$p_lrtnn[[i]] <- tout$p_value
  pardf$stat_lrtnn[[i]] <- tout$statistic
  pardf$df_lrtnn[[i]] <- tout$df

  ## old chisq ----
  nout <- chisq_gl4(gl = gl, g1 = g1, g2 = g2)
  pardf$p_chisq[[i]] <- nout$p_value
  pardf$df_chisq[[i]] <- nout$df
  pardf$stat_chisq[[i]] <- nout$statistic

  ## old polymapr ----
  gp <- exp(gl - apply(X = gl, MARGIN = 1, FUN = updog::log_sum_exp))
  pout <- polymapr_test(x = gp, g1 = g1, g2 = g2, type = "polymapR")
  pardf$p_polymapr[[i]] <- pout$p_value

  ## Bayes test ----
  trash <- capture.output(
    bout <- bayes_men_gl4(gl = gl, g1 = g1, g2 = g2, chains = 1)
  )
  pardf$lbf[[i]] <- bout$lbf
  pardf$alpha_bayes[[i]] <- bout$alpha
  pardf$xi1_bayes[[i]] <- bout$xi1
  pardf$xi2_bayes[[i]] <- bout$xi2
}

write_csv(x = pardf, file = "./output/sims/gl_altsims.csv")

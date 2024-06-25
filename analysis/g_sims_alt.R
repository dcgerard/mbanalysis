library(tidyverse)
library(menbayes)

pardf <- expand_grid(
  seed = 1:200,
  n = c(20, 200))

qvec <- rep(1/5, 5)

## new lrt params
pardf$p_lrt <- NA_real_
pardf$stat_lrt <- NA_real_
pardf$df_lrt <- NA_real_
pardf$alpha_lrt <- NA_real_
pardf$xi1_lrt <- NA_real_
pardf$xi2_lrt <- NA_real_

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

## other bayes priors
pardf$lbf_1 <- NA_real_
pardf$lbf_2 <- NA_real_
pardf$lbf_3 <- NA_real_
pardf$lbf_4 <- NA_real_

for (i in seq_len(nrow(pardf))) {
  cat(i, " of ", nrow(pardf), "\n")
  g1 <- 2
  g2 <- 2
  x <- c(stats::rmultinom(n = 1, size = pardf$n[[i]], prob = qvec))

  ## new lrt ----
  tout <- lrt_men_g4(x = x, g1 = g1, g2 = g2)
  pardf$p_lrt[[i]] <- tout$p_value
  pardf$stat_lrt[[i]] <- tout$statistic
  pardf$df_lrt[[i]] <- tout$df
  pardf$alpha_lrt[[i]] <- tout$alpha
  pardf$xi1_lrt[[i]] <- tout$xi1
  pardf$xi2_lrt[[i]] <- tout$xi2

  ## old chisq ----
  nout <- chisq_g4(x = x, g1 = g1, g2 = g2)
  pardf$p_chisq[[i]] <- nout$p_value
  pardf$df_chisq[[i]] <- nout$df
  pardf$stat_chisq[[i]] <- nout$statistic

  ## old polymapr ----
  pout <- polymapr_test(x = x, g1 = g1, g2 = g2, type = "polymapR")
  pardf$p_polymapr[[i]] <- pout$p_value

  ## new bayes ----
  trash <- capture.output(
    bout <- bayes_men_g4(x = x, g1 = g1, g2 = g2, chains = 1)
  )
  pardf$lbf[[i]] <- bout$lbf
  pardf$alpha_bayes[[i]] <- bout$alpha
  pardf$xi1_bayes[[i]] <- bout$xi1
  pardf$xi2_bayes[[i]] <- bout$xi2

  # Other Prior Settings ----
  trash <- capture.output(
    bout_1 <- bayes_men_g4(
      x = x,
      g1 = g1,
      g2 = g2,
      chains = 1,
      ts1 = 1/2,
      ts2 = 1/2,
      shape1 = 1,
      shape2 = 2)
  )
  pardf$lbf_1[[i]] <- bout_1$lbf

  trash <- capture.output(
    bout_2 <- bayes_men_g4(
      x = x,
      g1 = g1,
      g2 = g2,
      chains = 1,
      ts1 = 1/2,
      ts2 = 1/2,
      shape1 = 1/3,
      shape2 = 2/3)
  )
  pardf$lbf_2[[i]] <- bout_2$lbf

  trash <- capture.output(
    bout_3 <- bayes_men_g4(
      x = x,
      g1 = g1,
      g2 = g2,
      chains = 1,
      ts1 = 2,
      ts2 = 2,
      shape1 = 1,
      shape2 = 2)
  )
  pardf$lbf_3[[i]] <- bout_3$lbf

  trash <- capture.output(
    bout_4 <- bayes_men_g4(
      x = x,
      g1 = g1,
      g2 = g2,
      chains = 1,
      ts1 = 2,
      ts2 = 2,
      shape1 = 1/3,
      shape2 = 2/3)
  )
  pardf$lbf_4[[i]] <- bout_4$lbf

}

write_csv(x = pardf, file = "./output/sims/g_altsims.csv")

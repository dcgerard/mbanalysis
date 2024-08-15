library(tidyverse)
library(menbayes)
library(doFuture)
library(doRNG)

## Set up parallelization ----
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

# Create data frame of alternatives
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
df_alt$name <- paste0("alt_", 1:nrow(df_alt))

## Create simulation data frame
pardf <- expand_grid(
  seed = 1:200,
  n = c(20, 200),
  alt = c("unif", "random", df_alt$name),
  rd = 10)

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

## other bayes priors
pardf$lbf_1 <- NA_real_
pardf$lbf_2 <- NA_real_
pardf$lbf_3 <- NA_real_
pardf$lbf_4 <- NA_real_

chisq_gl_unknown_parents <- function(gl) {
  pval<- 0
  lfinal <- list(p_value = 0, g1 = NA, g2 = NA, df = NA, statistic = NA)
  for (g1 in 0:4) {
    for (g2 in g1:4) {
      suppressWarnings(
        lnow <- chisq_gl4(gl = gl, g1 = g1, g2 = g2)
      )
      if (lnow$p_value >= pval) {
        lfinal <- lnow
        lfinal$g1 <- g1
        lfinal$g2 <- g2
        pval <- lfinal$p_value
      }
    }
  }
  return(lfinal)
}

outdf <- foreach(
  i = seq_len(nrow(pardf)),
  .combine = rbind,
  .export = c("pardf")) %dorng% {
  set.seed(seed = pardf$seed[[i]])

  if (pardf$alt[[i]] == "unif") {
    qvec <- rep(1/5, 5)
  } else if (pardf$alt[[i]] == "random") {
    qvec <- stats::rexp(n = 5, rate = 1)
    qvec <- qvec / sum(qvec)
  } else {
    p0 <- df_alt$p0[df_alt$name == pardf$alt[[i]]]
    p1 <- df_alt$p1[df_alt$name == pardf$alt[[i]]]
    p2 <- df_alt$p2[df_alt$name == pardf$alt[[i]]]
    pvec <- c(p0, p1, p2)
    qvec <- stats::convolve(pvec, rev(pvec), type = "open")
    qvec[qvec < 0] <- 0 # for -1e-17
  }

  x <- c(stats::rmultinom(n = 1, size = pardf$n[[i]], prob = qvec))
  gl <- hwep::simgl(nvec = x, rdepth = pardf$rd[[i]], ret = "gl")

  ## new lrt ----
  tout <- lrt_men_gl4(gl = gl, g1 = NULL, g2 = NULL)
  g1 <- tout$p1
  g2 <- tout$p2
  pardf$p_lrt[[i]] <- tout$p_value
  pardf$stat_lrt[[i]] <- tout$statistic
  pardf$df_lrt[[i]] <- tout$df
  pardf$alpha_lrt[[i]] <- tout$alpha
  pardf$xi1_lrt[[i]] <- tout$xi1
  pardf$xi2_lrt[[i]] <- tout$xi2

  ## Basic LRT ----
  tout <- lrt_men_gl4(gl = gl, g1 = NULL, g2 = NULL, pp = FALSE, dr = FALSE, alpha = 0, xi1 = 1/3, xi2 = 1/3)
  pardf$p_lrtnn[[i]] <- tout$p_value
  pardf$stat_lrtnn[[i]] <- tout$statistic
  pardf$df_lrtnn[[i]] <- tout$df

  ## old chisq ----
  nout <- chisq_gl_unknown_parents(gl = gl)
  pardf$p_chisq[[i]] <- nout$p_value
  pardf$df_chisq[[i]] <- nout$df
  pardf$stat_chisq[[i]] <- nout$statistic

  ## old polymapr ----
  gp <- exp(gl - apply(X = gl, MARGIN = 1, FUN = updog::log_sum_exp))
  pout <- polymapr_test(x = gp, g1 = NULL, g2 = NULL, type = "menbayes")
  pardf$p_polymapr[[i]] <- pout$p_value

  ## Bayes test ----
  trash <- capture.output(
    bout <- bayes_men_gl4(gl = gl, g1 = g1, g2 = g2, chains = 1)
  )
  pardf$lbf[[i]] <- bout$lbf
  pardf$alpha_bayes[[i]] <- bout$alpha
  pardf$xi1_bayes[[i]] <- bout$xi1
  pardf$xi2_bayes[[i]] <- bout$xi2

  ## Other Prior Settings ----
  trash <- capture.output(
    bout_1 <- bayes_men_gl4(
      gl = gl,
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
    bout_2 <- bayes_men_gl4(
      gl = gl,
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
    bout_3 <- bayes_men_gl4(
      gl = gl,
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
    bout_4 <- bayes_men_gl4(
      gl = gl,
      g1 = g1,
      g2 = g2,
      chains = 1,
      ts1 = 2,
      ts2 = 2,
      shape1 = 1/3,
      shape2 = 2/3)
  )
  pardf$lbf_4[[i]] <- bout_4$lbf

  pardf[i, ]
}

## Unregister workers ----
if (nc > 1) {
  plan(sequential)
}

write_csv(x = outdf, file = "./output/sims/gl_altsims.csv")

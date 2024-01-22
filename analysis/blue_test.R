library(updog)
library(menbayes)
library(tidyverse)
library(future)
library(foreach)
library(doFuture)
library(rngtools)
library(doRNG)
bluefits <- readRDS("./output/blue/bluefits.RDS")

## Determine number of cores ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nc <- 1
} else {
  eval(parse(text = args[[1]]))
}
cat(nc, "\n")

## Register workers ----
if (nc == 1) {
  registerDoSEQ()
  plan(sequential)
} else {
  registerDoFuture()
  registerDoRNG()
  plan(multisession, workers = nc)
  if (getDoParWorkers() == 1) {
    stop("nc > 1, but only one core registered")
  }
}

## This will give you an array with dimensions SNPs by Individuals by Genotypes
gl <- format_multidog(bluefits, varname = paste0("logL_", 0:4))
p1mat <- gl[, "indigocrisp", ]
p2mat <- gl[, "sweetcrisp", ]
gl <- gl[, setdiff(dimnames(gl)[[2]], c("indigocrisp", "sweetcrisp")), ]

##Build dataframe
blue_df <- data.frame(snp = dimnames(gl)[[1]])
blue_df$p1 <- NA_real_
blue_df$p2 <- NA_real_
blue_df$lbf <- NA_real_
blue_df$pm_alpha <- NA_real_
blue_df$pm_xi1 <- NA_real_
blue_df$pm_xi2 <- NA_real_
blue_df$chisq_stat <- NA_real_
blue_df$chisq_pvalue <- NA_real_
blue_df$lrt_df <- NA_real_
blue_df$lrt_stat <- NA_real_
blue_df$lrt_pvalue <- NA_real_
blue_df$lrt_alpha <- NA_real_
blue_df$lrt_xi1 <- NA_real_
blue_df$lrt_xi2 <- NA_real_
blue_df$lrtnn_stat <- NA_real_
blue_df$lrtnn_pvalue <- NA_real_
blue_df$polymapr_pvalue <- NA_real_

## Sanity checks
stopifnot(nrow(blue_df) == dim(gl)[[1]])

outdf <- foreach(i = seq_len(dim(gl)[[1]]), .combine = rbind, .export = c("blue_df")) %dopar% {
  glmat <- gl[i , ,]
  goodind <- apply(glmat, 1, function(x) all(!is.na(x)))
  glmat <- glmat[goodind, , drop = FALSE]
  p1 <- p1mat[i, ]
  p2 <- p2mat[i, ]

  ## Fit Bayes test here
  bout <- bayes_men_gl4(gl = glmat, g1 = p1, g2 = p2, chains = 1)

  blue_df$lbf <- bout$lbf
  blue_df$pm_alpha <- bout$alpha
  blue_df$pm_xi1 <- bout$xi1
  blue_df$pm_xi2 <- bout$xi2

  ## Fit LRT
  lout <- lrt_men_gl4(gl = glmat, g1 = p1, g2 = p2)
  blue_df$p1 <- lout$p1
  blue_df$p2 <- lout$p2
  blue_df$lrt_stat <- lout$statistic
  blue_df$lrt_df <- lout$df
  blue_df$lrt_pvalue <- lout$p_value
  blue_df$lrt_alpha <- lout$alpha
  blue_df$lrt_xi1 <- lout$xi1
  blue_df$lrt_xi2 <- lout$xi2

  ## Fit LRT nn
  lnnout <- lrt_men_gl4(gl = glmat, g1 = p1, g2 = p2, pp = FALSE, dr = FALSE)
  blue_df$lrtnn_stat <- lnnout$statistic
  blue_df$lrtnn_pvalue <- lnnout$p_value

  ## Fit chi-squared test
  if (!(lout$p1 %in% c(0, 4) && lout$p2 %in% c(0, 4))) {
    cout <- chisq_gl4(gl = glmat, g1 = lout$p1, g2 = lout$p2)
    blue_df$chisq_pvalue <- cout$p_value
    blue_df$chisq_stat <- cout$statistic
  }

  ## Fit polmapr test
  gp <- exp(glmat - apply(X = glmat, MARGIN = 1, FUN = updog::log_sum_exp))
  pout <- polymapr_test(x = gp, g1 = lout$p1, g2 = lout$p2, type = "menbayes")
  blue_df$polymapr_pvalue[[i]] <- pout$p_value

  blue_df[i, ]
}

## Unregister workers ----
if (nc > 1) {
  plan(sequential)
}

write.csv(outdf, "./output/blue/blue_df.csv")

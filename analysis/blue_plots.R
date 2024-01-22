library(tidyverse)
library(updog)
library(menbayes)
library(GGally)
library(gridExtra)
library(ggpubr)
library(xtable)

bdf <- read_csv("./output/blue/blue_df.csv")
uout <- readRDS("./output/blue/bluefits.RDS")

## Add parent info for plotting purposes
uout$inddf |>
  select(snp, ind, p1ref = ref, p1size = size) |>
  filter(ind == "sweetcrisp") |>
  select(-ind) ->
  df_sc

uout$inddf |>
  select(snp, ind, p2ref = ref, p2size = size) |>
  filter(ind == "indigocrisp")  |>
  select(-ind) ->
  df_ic

uout$snpdf$p1ref <- NULL
uout$snpdf$p1size <- NULL
uout$snpdf$p2ref <- NULL
uout$snpdf$p2size <- NULL

uout$snpdf |>
  left_join(df_sc, by = join_by(snp)) |>
  left_join(df_ic, by = join_by(snp)) ->
  uout$snpdf

## Filter monomorphic SNPs ----
mdf <- tibble(snp = uout$snpdf$snp, mono = ismono <- apply(uout$snpdf[paste0("Pr_", 0:4)], 1, max) > 0.95)
mdf |>
  filter(!mono) ->
  mdf

## Number of non-monomorphic SNPs ----
sum(!mdf$mono)

## Keep only main linkage groups ----
bdf |>
  semi_join(mdf, by = join_by(snp)) |>
  filter(!(p1 %in% c(0, 4) & p2 %in% c(0, 4))) |>
  mutate(Chromosome = parse_number(str_extract(snp, "^\\d+\\_")),
         Position = parse_number(str_extract(snp, "\\_\\d+$"))) |>
  filter(Chromosome <= 22) ->
  bdf_nm

## Number of remaining SNPs ----
table(bdf_nm$Chromosome)
nrow(bdf_nm)

## Number of "bad" SNPs ----
mean(bdf_nm$lbf < -16)
mean(p.adjust(p = bdf_nm$lrt_pvalue, method = "bonferroni") < 0.05)
mean(p.adjust(p = bdf_nm$polymapr_pvalue, method = "bonferroni") < 0.05)
mean(p.adjust(p = bdf_nm$chisq_pvalue, method = "bonferroni") < 0.05, na.rm = TRUE)

## Plot of p-values and bayes Factors ----
bdf_nm |>
  select(`log-BF` = lbf, `LRT log p-value` = lrt_pvalue, `polymapR log p-value` = polymapr_pvalue) |>
  mutate(`LRT log p-value` = log(`LRT log p-value`), `polymapR log p-value` = log(`polymapR log p-value`)) |>
  ggpairs(upper = list(continuous = "blank"),
    lower = list(continuous = "points"),
    diag = list(continuous = "blankDiag")) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) ->
  pl

#' Function from StackOverflow
#' https://stackoverflow.com/questions/42654928/how-to-show-only-the-lower-triangle-in-ggpairs
gpairs_lower <- function(g){
  g$plots <- g$plots[-(1:g$nrow)]
  g$yAxisLabels <- g$yAxisLabels[-1]
  g$nrow <- g$nrow -1

  g$plots <- g$plots[-(seq(g$ncol, length(g$plots), by = g$ncol))]
  g$xAxisLabels <- g$xAxisLabels[-g$ncol]
  g$ncol <- g$ncol - 1

  g
}

pdf(file = "./output/blue/plots/pairs.pdf", height = 4, width = 4, family = "Times")
gpairs_lower(pl)
dev.off()

## Plot SNPs where LRT indicates reject but polymapR indicates no reject ----
bdf_nm |>
  filter(lrt_pvalue > 0.99, polymapr_pvalue < 0.01) |>
  arrange(polymapr_pvalue) |>
  select(snp, p1, p2, lbf, lrt_pvalue, polymapr_pvalue) ->
  polystrong

usub <- filter_snp(uout, snp %in% polystrong$snp[1:5])
plist <- plot(usub, indices = 1:5)
for (i in seq_along(plist)) {
  plist[[i]] <- plist[[i]] +
    scale_alpha(range = c(1, 1), guide = "none")
  if (i == 1) {
    leg <- get_legend(plist[[1]])
  }
  plist[[i]] <- plist[[i]] +
    theme(legend.position = "none")
}
plist[[6]] <- leg
ml <- marrangeGrob(plist, nrow = 3, ncol = 2, top = NULL)
ggsave(
  filename = "./output/blue/plots/polystrong.pdf",
  plot = ml,
  width = 4,
  height = 6,
  family = "Times")


## lrt reject, polymapr no reject ----
bdf_nm |>
  filter(lrt_pvalue < 0.05, polymapr_pvalue > 0.95) |>
  arrange(lrt_pvalue) |>
  select(snp, p1, p2, lbf, lrt_pvalue, polymapr_pvalue)  ->
  lrtstrong

usub <- filter_snp(uout, snp %in% lrtstrong$snp[1:5])
plist <- plot(usub, indices = 1:5)
for (i in seq_along(plist)) {
  plist[[i]] <- plist[[i]] +
    scale_alpha(range = c(1, 1), guide = "none")
  if (i == 1) {
    leg <- get_legend(plist[[1]])
  }
  plist[[i]] <- plist[[i]] +
    theme(legend.position = "none")
}
plist[[6]] <- leg
ml <- marrangeGrob(plist, nrow = 3, ncol = 2, top = NULL)
ggsave(
  filename = "./output/blue/plots/lrtstrong.pdf",
  plot = ml,
  width = 4,
  height = 6,
  family = "Times")

## Table of polystrong and lrtstrong top SNPs
bind_rows(lrtstrong[1:5, ], polystrong[1:5, ]) |>
  rename(LRT = lrt_pvalue, polymapR = polymapr_pvalue, Bayes = lbf, SNP = snp, ell1 = p1, ell2 = p2) |>
  xtable(
    display = rep("g", ncol(lrtstrong) + 1),
    label = "tab:blue.diff",
    caption = 'In the first five SNPs (plotted in Figure \\ref{fig:lrt.strong}), the LRT indicates segregation distortion while the \\textsf{polymapR} test indicates no segregation distortion, while in the last five SNPs (plotted in Figure \\ref{fig:poly.strong}) the LRT indicates no segregation distortion while the \\textsf{polymapR} test indicates segregation distortion. Parent genotypes (ell1 and ell2) are listed, along with the log Bayes factors (Bayes), the LRT $p$-values ("LRT"), and the \\textsf{polyampR} $p$-values ("polymapR").') |>
  print(include.rownames = FALSE,
        file = "./output/blue/plots/tensnps.tex")

## re-check those SNPs with polymapR
gp <- format_multidog(x = uout, varname = paste0("Pr_", 0:4))
gp <- gp[, setdiff(dimnames(gp)[[2]], c("sweetcrisp", "indigocrisp")), ]
gl <- format_multidog(x = uout, varname = paste0("logL_", 0:4))
gl <- gl[, setdiff(dimnames(gl)[[2]], c("sweetcrisp", "indigocrisp")), ]
genomat <- format_multidog(x = uout, varname = "geno")

i <- 5
gpmat <- gp[polystrong$snp[[i]], , ]
glmat <- gl[polystrong$snp[[i]], , ]
geno <- unlist(genomat[polystrong$snp[[i]], ])
polymapr_test(x = gpmat, g1 = polystrong$p1[[i]], g2 = polystrong$p2[[i]], type = "polymapR")$p_value
polymapr_test(x = gpmat, g1 = polystrong$p1[[i]], g2 = polystrong$p2[[i]], type = "menbayes")$p_value

lrt_men_gl4(gl = glmat, g1 = polystrong$p1[[i]], g2 = polystrong$p2[[i]])$p_value

table(geno)

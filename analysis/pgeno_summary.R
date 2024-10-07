library(updog)
bout <- readRDS(file = "./output/blue/bluefits.RDS")
geno <- format_multidog(bout, "geno")
gsub <- geno[, colnames(geno) %in% c("indigocrisp", "sweetcrisp")]

gsub[gsub == 3] <- 1
gsub[gsub == 4] <- 0

tab <- table(gsub[, "indigocrisp"], gsub[, "sweetcrisp"])
dimnames(tab) <- list("Parent1" = c("nullplex", "simplex", "duplex"),
                      "Parent2" = c("nullplex", "simplex", "duplex"))
tab

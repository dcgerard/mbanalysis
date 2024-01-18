#################
## Blueberry Analysis
#################
library(updog)
library(future)
load("./data/updog_input_240ind_Sweet_Indi.Rdata")

## Get number of cores from makefile
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nc <- 1
} else {
  eval(parse(text = args[[1]]))
}
cat(nc, "\n")

## Verify that cbind won't frell us
stopifnot(rownames(ao.sibling) == names(ao.sweetcrisp),
          rownames(ao.sibling) == names(ao.indigocrisp))
stopifnot(rownames(dp.sibling) == names(dp.sweetcrisp),
          rownames(dp.sibling) == names(dp.indigocrisp))

altmat <- cbind(ao.sweetcrisp, ao.indigocrisp, ao.sibling)
colnames(altmat)[1:2] <- c("sweetcrisp", "indigocrisp")
altmat <- as.matrix(altmat)
sizemat <- cbind(dp.sweetcrisp, dp.indigocrisp, dp.sibling)
colnames(sizemat)[1:2] <- c("sweetcrisp", "indigocrisp")
sizemat <- as.matrix(sizemat)

## Run multidog
mout <- multidog(refmat = altmat,
                 sizemat = sizemat,
                 ploidy = 4,
                 model = "f1pp",
                 p1_id = "sweetcrisp",
                 p2_id = "indigocrisp",
                 nc = nc)

saveRDS(object = mout, file = "./output/blue/bluefits.RDS")

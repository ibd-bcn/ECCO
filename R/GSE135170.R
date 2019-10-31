library("GEOquery")
gse <- read.table("data/GSE135170_colonoid.072919.txt")

meta <- strsplit(colnames(gse), split = "_|\\.")
meta <- t(simplify2array(meta))
meta <- meta[, 2:3]
colnames(meta) <- c("Origin", "Type")

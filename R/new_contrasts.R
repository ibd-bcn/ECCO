# Loading libraries ####
library("readxl")
library("edgeR")
library("limma")
library("dplyr")

# Meta data ####
folder <- "data"

sam_xlsx <- file.path(folder, "Final_samples.xlsx")
sam <- read_xlsx(path = sam_xlsx, range = "B1:X162")

contr_xlsx <- file.path(folder, "contrast_characteristics_original.xlsx")
contr <- read_xlsx(path = contr_xlsx, range = "H1:J37")
colnames(contr)[1] <- "contrast"

contr <- contr %>%
    mutate(contrast = gsub("/", " ", contrast),
           `0` = gsub("/", " ", `0`))

meta <- sam %>%
    mutate(reanalyzed = endsWith(`AGILENT/ARRAY CODE`, "V"),
           colname = `AGILENT/ARRAY CODE`,
           development = tolower(`STEM/DIFF`),
           status = `SAMPLE STATUS`)

# RNAseq preprocessing ####
counts <- read.table(here::here("data", "ORGANOIDS_STAR_RSEM_GENES.txt"),
                     check.names = FALSE)

samples <- gsub("_.*", "", colnames(counts))
uniq_samples <- unique(samples)
genes <- gsub("\\..*", "", rownames(counts))

# To summarize the samples together
group_samples <- function(expr, names) {
    umat <- matrix(nrow = nrow(expr), ncol = length(names))
    for (i in 1:length(names)) {
        umat[, i] <- rowSums(expr[, which(samples == names[i])])
    }
    colnames(umat) <- names
    rownames(umat) <- rownames(expr)
    umat
}
umat <- group_samples(counts, uniq_samples)
umat <- umat[rowSums(umat) != 0, ]

# Subset to just old samples
expr <- umat[, colnames(umat) %in% meta$colname[meta$reanalyzed]]
meta <- meta[meta$colname %in% colnames(expr), ]
expr <- expr[, meta$`AGILENT/ARRAY CODE`]

# Creating contrasts ####
meta2 <- meta[, c("colname", "development", "GROUP", "TYPE", "LOCATION",
                  "status", "INVOLVEMENT")]
colnames(meta2) <- tolower(colnames(meta2))

comp <- lapply(contr$contrast, function(x, m){
    g1 <- unlist(strsplit(contr$`0`[contr$contrast == x], split = "_"),
           recursive = FALSE, use.names = FALSE)
    g2 <- unlist(strsplit(contr$`1`[contr$contrast == x], split = "_"),
           recursive = FALSE, use.names = FALSE)
    y <- vector("numeric", nrow(m))
    m2 <- matrix(m %in% tolower(g1), ncol = ncol(m), nrow = nrow(m))
    m3 <- matrix(m %in% tolower(g2), ncol = ncol(m), nrow = nrow(m))
    y[rowSums(m2) == length(g1)] <- 2
    y[rowSums(m3) == length(g2)] <- 1
    y
}, m = tolower(as.matrix(meta2[, 2:7])))
# * Formatting as table ####
s <- simplify2array(comp)
colnames(s) <- contr$contrast
rownames(s) <- meta$colname
s[s == 0] <- NA

write.csv(s, "processed/new_comparisons.csv")
# * Checking contrasts ####
con <- apply(s, 2, table)
# To check the contrasts
k <- vapply(con, function(x){any(x[c("1", "2")] == 1)}, logical(1L))
stopifnot(sum(lengths(con) < 2 | k | is.na(k)) <= 10)
incomplete <- lengths(con) < 3
con[incomplete] <- lapply(con[incomplete],
                          function(x){
                              diff_names <- c("1", "2")
                              y <- setdiff(diff_names, names(x))
                              z <- 0
                              names(z) <- y
                              z <- c(x, z)
                              z[diff_names]
                          })
values <- t(simplify2array(con))
write.csv(values, "processed/new_contrasts.csv")

# GETS ####
# perl ~/Documents/projects/GETS/gets.pl \
#  --matrix=./processed/old.tsv \
#  --geneinfo=./processed/genes_juanjo.tsv \
#  --sampleinfo=./processed/samples_info_old.tsv \
#  --colors=./processed/colors_info_old.tsv \
#  --output=./processed/GETS_old_juanjo_v5 \
#  --center=TRUE --overwrite=TRUE

library("edgeR")
library("tidyr")
library("ggplot2")
library("ggrepel")
library("dplyr")
library("sva")
library("here")
library("matrixStats")
library("hgu219.db") # Newest
library("hgu219hsentrezg.db") # From local storage...

## Metadata ####
rename_cols <- c(colname = "AGILENT/ARRAY CODE",
                 cell_type = "STEM/DIFF",
                 status = "SAMPLE STATUS",
                 AaD = "Age at diagnosis",
                 GROUP = "GROUP",
                 age = "age",
                 TYPE = "TYPE",
                 LOCATION = "LOCATION",
                 SAMPLE = "SAMPLE"
)

meta <- readxl::read_excel(here::here("data", "variables CD cohort.xlsx")) %>%
    dplyr::select(rename_cols) %>%
    filter(GROUP != "UC") %>%
    mutate(GROUP = toupper(GROUP),
           reanalyzed = endsWith(colname, "V"),
           # Age at diagnosis correted as per Isa instructions
           age = ifelse(colname == "34V", 52, age),
           # Isa: consider as adult above the 17 years
           TYPE = case_when(age <= 17 & !is.na(age) ~ "pediatric",
                            age > 17 & !is.na(age) ~ "adult",
                            TRUE ~ TYPE),
           SAMPLE2 = gsub("/.*", "", gsub(" .*", "", trimws(gsub("STEM|DIFF.?", "", gsub("RNA ", " ", SAMPLE)))))
    )

# RNAseq ####
counts <- read.table(here::here("data", "ORGANOIDS_STAR_RSEM_GENES.txt"), check.names = FALSE)

samples <- gsub("_.*", "", colnames(counts))
uniq_samples <- unique(samples)
genes <- gsub("\\..*", "", rownames(counts))

meta <- subset(meta, colname %in% uniq_samples)

# To summarize the samples together
group_samples <- function(expr, names) {
    umat <- matrix(nrow = nrow(expr), ncol = length(names))
    for (i in 1:length(names)) {
        umat[, i] <- apply(expr[, which(samples == names[i])], 1, sum)
    }
    colnames(umat) <- names
    rownames(umat) <- rownames(expr)
    umat
}
umat <- group_samples(counts, uniq_samples)
umat <- umat[rowSums(umat) != 0, ]

expr <- umat[, meta$colname[meta$reanalyzed]]

dge_rna <- DGEList(expr)
dge_rna <- calcNormFactors(dge_rna)
voom_rna <- voom(dge_rna, plot = FALSE, normalize.method = "quantile")

# microarray ####
microarrays <- read.table(
    here::here("data", "nov072019_Fulldata_organoides_nonFiltered.txt"),
    check.names = FALSE, row.names = 1, sep = "\t")
microarrays <- as.matrix(microarrays)
microarray_n <- read.csv(
    here::here("data", "nov072019_Fulldata_annotation.txt"),
    check.names = FALSE, row.names = 1, sep = "\t")

# Use the default normalization
# dge_micro <- DGEList(microarrays)
# dge_micro <- calcNormFactors(dge_micro)
# voom_micro <- voom(dge_micro, plot = FALSE, normalize.method = "quantile")

# Prepare the data ####
# Subset samples and reorder them
shared_names <- intersect(meta$SAMPLE, colnames(microarrays))
expr_rna <- voom_rna$E[, match(shared_names, meta$SAMPLE[meta$reanalyzed])]
expr_micro <- microarrays[, shared_names]
# expr_micro <- voom_micro$E[, shared_names]

# Subset genes
genes_rna <- data.frame(real = rownames(expr_rna),
                        hugo = gsub(".*_(.*)", "\\1", rownames(expr_rna)))

microarray_n <- microarray_n[microarray_n$probes %in% rownames(expr_micro), ]

shared_genes <- intersect(genes_rna$hugo, microarray_n$hugo)
shared_genes <- shared_genes[!is.na(shared_genes)]

microarray_n <- microarray_n[microarray_n$hugo %in% shared_genes, ]
genes_rna <- genes_rna[genes_rna$hugo %in% shared_genes, ]

microarray_n <- microarray_n[match(shared_genes, microarray_n$hugo), ]
microarray_n <- droplevels(microarray_n)
genes_rna <- genes_rna[match(shared_genes, genes_rna$hugo), ]

expr_rna <- expr_rna[genes_rna$real, ]
expr_micro <- expr_micro[microarray_n$probes, ]

colnames(expr_micro) <- colnames(expr_rna)
rownames(expr_rna) <- genes_rna$hugo
rownames(expr_micro) <- microarray_n$hugo

write.csv(expr_rna, file = "processed/RNAseq_samples.csv")
write.csv(expr_micro, file = "processed/microarray_samples.csv")

# Scale/normalize
r_values <- rowMedians(scale(expr_rna, center = TRUE, scale = FALSE))
m_values <- rowMedians(scale(expr_micro, center = TRUE, scale = FALSE))

names(r_values) <- rownames(expr_rna)
names(m_values) <- rownames(expr_micro)

# Compare ####
pars <- par(pty = "s")
plot(r_values, m_values, xlab = "RNA-seq", ylab = "Microarray",
     pch = 16, cex = 0.5, asp = 1, main = "Comparing organoids expression")
abline(a = 0, b = 1, col = "red", lty = 2)
# text(x = 1, y = -1, labels = "1446 genes")
par(pars)

data.frame(rnaseq = r_values,
           micro = m_values) %>%
    ggplot(aes(rnaseq, micro)) +
    geom_point() +
    geom_smooth(method = "lm") +
    labs(x = "RNA-seq", y = "Microarray",
         title = "Comparing organoids expression") +
    coord_fixed(xlim = c(-7, 7), ylim = c(-7, 7)) +
    theme_minimal()
cor.test(r_values, m_values)
lm(r_values ~ m_values) %>%
    broom::tidy()

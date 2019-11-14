library("affy")
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


# Functions  ####
# taken from correlation_CEL_vs_RNAseq.R
remove.genes.with.low.expression <- function(data, q) {
  data.means <- as.data.frame(rowMeans(data)) # RNA-seq data: removing genes with low expression ( < 10th quantile)
  min.exp <- quantile(data.means[, 1], q, na.rm = TRUE)
  # data <- data[rownames(data.means)[data.means[,1] > min.exp],]
  data <- data[data.means[, 1] > min.exp, ]

  if (nrow(data) == 0) {
    stop("[     ERROR]\tNo genes in data\n")
  }

  return(data)
}

remove.genes.with.low.CV <- function(data, q) {
  cv <- apply(data, 1, sd) / apply(data, 1, mean)
  cut <- quantile(cv, q)
  filter <- which(cv > cut)

  data <- data[filter, ]

  if (nrow(data) == 0) {
    stop("[     ERROR]\tNo genes in data\n")
  }

  return(data)
}

plot.CEL.vs.RNAseq <- function(CEL.data, RNAseq.data, file, xlab, ylab) {
  cel <- rowMeans(CEL.data)
  rnaseq <- rowMeans(RNAseq.data)

  min <- min(cel, rnaseq)
  max <- max(cel, rnaseq)

  gene.list <- intersect(rownames(CEL.data), rownames(RNAseq.data))

  pdf(file)

  x.list <- vector(, length(gene.list))
  y.list <- vector(, length(gene.list))

  plot(1, type = "n", xlim = c(min, max), ylim = c(min, max), main = "CEL vs RNA-seq\nComparison of means", xlab = xlab, ylab = ylab)

  for (i in 1:length(gene.list))
  {
    gene <- gene.list[i]

    x <- mean(RNAseq.data[gene, ])
    y <- mean(CEL.data[gene, ])

    points(x, y, pch = 20, cex = 0.6, col = "blue")

    x.list[i] <- x
    y.list[i] <- y
  }

  segments(min, min, max, max, lty = "dashed")

  c <- cor.test(x.list, y.list, method = "pearson")

  legend("topright", legend = paste("r =", round(c$estimate, 2)), cex = 1, bty = "n")

  dev.off()
}

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

# microarray ####
# Already normalized by RMA
microarrays <- readRDS(here::here("data", "microarrays.RDS"))
microarray_n <- read.csv(here::here("data", "nov142019_Fulldata_annotation.txt"),
                         check.names = FALSE, row.names = 1, sep = "\t")

shared_names <- intersect(meta$SAMPLE, colnames(microarrays))
expr <- expr[, meta$colname[meta$SAMPLE %in% shared_names]]
microarrays <- microarrays[, shared_names]
micro <- affy::rma(microarrays, verbose = FALSE, normalize = TRUE)


dge_rna <- DGEList(expr)
dge_rna <- calcNormFactors(dge_rna)
voom_rna <- voom(dge_rna, plot = FALSE, normalize.method = "quantile")
voom_rna <- voom_rna$E

# Use the default normalization ####
low_CV <- 0.25
low_expression <- 0.1
e_micro <- remove.genes.with.low.expression(exprs(micro), low_expression)
expr_micro <- remove.genes.with.low.CV(e_micro, low_CV)


e_rna <- remove.genes.with.low.expression(voom_rna, low_expression)
expr_rna <- remove.genes.with.low.CV(e_rna, low_CV)

# Scale and do the median
r_values <- scale(rowMedians(expr_rna), center = TRUE, scale = FALSE)
m_values <- scale(rowMedians(expr_micro), center = TRUE, scale = FALSE)

dim(r_values) <- NULL
dim(m_values) <- NULL

names(r_values) <- rownames(expr_rna)
names(m_values) <- rownames(expr_micro)

# Prepare the data ####
# Subset genes
genes_rna <- data.frame(real = rownames(expr_rna),
                        hugo = gsub(".*_(.*)", "\\1", rownames(expr_rna)))

microarray_n <- microarray_n[microarray_n$probes %in% rownames(expr_micro), ]

shared_genes <- intersect(genes_rna$hugo, microarray_n$hugo)
stopifnot(length(shared_genes) > 1)
shared_genes <- shared_genes[!is.na(shared_genes)]

microarray_n <- microarray_n[microarray_n$hugo %in% shared_genes, ]
genes_rna <- genes_rna[genes_rna$hugo %in% shared_genes, ]

microarray_n <- microarray_n[match(shared_genes, microarray_n$hugo), ]
microarray_n <- droplevels(microarray_n)
genes_rna <- genes_rna[match(shared_genes, genes_rna$hugo), ]

r_values <- r_values[genes_rna$real]
m_values <- m_values[microarray_n$probes]

# Compare ####
cors <- cor.test(r_values, m_values, method = "spearman")
# L'escalat sembla com si es fes la mitjana per mostra no per gen!!!
# Llavors la correlació s'aproxima (tot i que els números no)
pars <- par(pty = "s")
plot(r_values, m_values, xlab = "RNA-seq", ylab = "Microarray",
     pch = 16, cex = 0.5, asp = 1, main = "Comparing organoids expression")
legend("topleft", legend = paste("r", signif(cors$estimate, 1)))
abline(a = 0, b = 1, col = "red", lty = 2)
# text(x = 1, y = -1, labels = "1446 genes")
par(pars)

data.frame(rnaseq = r_values,
           micro = m_values) %>%
    ggplot(aes(rnaseq, micro)) +
    ggpointdensity::geom_pointdensity() +
    viridis::scale_color_viridis() +
    geom_smooth(method = "lm") +
    labs(x = "RNA-seq", y = "Microarray",
         title = "Comparing organoids expression") +
    coord_fixed(xlim = c(-7, 7), ylim = c(-7, 7)) +
    theme_minimal()
lm(r_values ~ m_values) %>%
    broom::tidy()

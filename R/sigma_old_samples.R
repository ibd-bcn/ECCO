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
library("sva")

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
expr <- expr[, meta$colname]
stopifnot(meta$colname == colnames(expr)) # Same order
sigma <- meta$LOCATION == "sigma"
meta <- meta[sigma, ]
expr <- expr[, sigma]

# Design ####
grouping <- apply(meta[, c("cell_type", "GROUP", "TYPE"), drop = FALSE], 1,
                  paste0, collapse = "_")
grouping <- as.factor(grouping)
design <- model.matrix(~0 + grouping, grouping)
colnames(design) <- levels(grouping)
rownames(design) <- meta$colname


# DIFF eq 1, STEM to -1
contrasts <- c(
    diff_stem = contrasting("^DIFF", "^STEM", design), # Siempre por separado
    CD_CTRL = contrasting("_CD_", "_CTRL_", design),
    pediatric_adult = contrasting("pediatric", "adult", design)
)

contr.matrix <- makeContrasts(contrasts = contrasts, levels = colnames(design))
colnames(contr.matrix) <- names(contrasts)


pick <- design %*% contr.matrix # To know which samples are taken
# Look if there is any contrast with just one sample on one side

dge_rna <- DGEList(expr)
keep <- filterByExpr(dge_rna)
dge_rna <- dge_rna[keep, ]


dge_rna <- calcNormFactors(dge_rna, method = "TMM")
# Limma p 125-126 two vooms
voom_rna <- voom(dge_rna, plot = TRUE, normalize.method = "cyclicloess")
vfit <- lmFit(voom_rna, design = design)

# Continue evaluation whole ####
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)

l <- lapply(seq_along(contrasts), function(x) {
    topTable(efit, coef = x, number = Inf)
})

dt <- decideTests(efit, lfc = log2(1.5), adjust.method = "fdr")
dtt <- summary(dt)
dtt


# Prepare for GET ####
# perl ~/Documents/projects/GETS/gets.pl --matrix=./processed/old.tsv --geneinfo=./processed/gene_info_old.tsv --sampleinfo=./processed/STEM_samples_info_old.tsv --colors=./processed/colors_info_old.tsv --output=./processed/STEM_GETS_old --center=TRUE --overwrite=TRUE
genes <- gsub("\\..*", "", rownames(dt))
genes_n <- gsub(".*\\.[0-9]*_", "", rownames(dt))
gene_names <- mapIds(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", "SYMBOL")
# Possible step to filter the genes to reduce the number of rows
# summary(is.na(gene_info$gene_names))
# kruskal.test(list(table(is_na), table(is_pr))) Is not conclusive
gene_info <- data.frame(rownames = rownames(dt), genes, gene_names, genes_n)
gene_info$Description <- mapIds(org.Hs.eg.db,
                                keys = as.character(gene_info$genes),
                                keytype = "ENSEMBL", column = "GENENAME")

# Adding data from the comparisons

tt_all <- lapply(seq_along(contrasts), function(x) {
    tt <- topTable(efit, coef = x, number = Inf, adjust.method = "BH")
    sign <- sign(tt$logFC)
    tt$FC <- sign*2^(abs(tt$logFC))
    above <- abs(tt$FC) > 1.5
    signif <- tt$P.Value < 0.05
    signiff <- tt$adj.P.Val < 0.05
    tt$sign <- ""

    tt$sign[signif & sign < 0 & above] <- "DW"
    tt$sign[signiff & sign < 0 & above] <- "DDW"
    tt$sign[signif & sign > 0 & above] <- "UP"
    tt$sign[signiff & sign > 0 & above] <- "UUP"
    tt <- tt[, !grepl("B$|AveExpr$|t$", colnames(tt))]
    colnames(tt) <- paste0("contrast", x, "_", names(contrasts)[x], "_", colnames(tt))
    tt
})

tt_all_m <- do.call(cbind, tt_all)
sign_c <- grep("_sign$", colnames(tt_all_m))
FC_c <- grep("_FC$", colnames(tt_all_m))
adj_c <- grep("adj.P.Value$", colnames(tt_all_m))
raw_c <- grep("P.Value$", colnames(tt_all_m))

tt_all_m <- tt_all_m[, c(sign_c, FC_c, adj_c, raw_c)]
gene_info2 <- cbind(gene_info, tt_all_m)


write.table(gene_info2[, -1], "processed/gene_info_old_sigma.tsv", sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE, na = "")

stopifnot(ncol(voom_rna$E) == nrow(meta))
r <- order(meta$TYPE, meta$cell_type, meta$LOCATION, meta$GROUP)
meta_o <- meta[r, -10]
meta_o <- meta_o[, c("colname", "TYPE", "cell_type", "GROUP", "status",
                     "age", "AaD", "SAMPLE", "SAMPLE2")]
expr_csv <- voom_rna$E[, r]
stopifnot(ncol(expr_csv) == nrow(meta))
stopifnot(nrow(gene_info2) == nrow(expr_csv))
expr2_csv <- cbind.data.frame(genes = rownames(expr_csv), expr_csv)
write.table(expr2_csv, file = "processed/old_sigma.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
stem <- cbind.data.frame(genes = rownames(expr_csv), expr_csv[, meta_o$cell_type == "STEM"])
write.table(stem,
            file = "processed/STEM_old_sigma.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")
diff <- cbind.data.frame(genes = rownames(expr_csv), expr_csv[, meta_o$cell_type == "DIFF"])
write.table(diff,
            file = "processed/DIFF_old_sigma.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")

meta_o %>%
    dplyr::filter(colname %in% colnames(expr2_csv)) %>%
    dplyr::select(-colname) %>%
    t() %>%
    write.table(file = "processed/samples_info_old_sigma.tsv", sep = "\t",
                row.names = TRUE, quote = FALSE, col.names = FALSE, na = "")
meta_o %>%
    dplyr::filter(colname %in% colnames(stem)[-1]) %>%
    dplyr::select(-colname) %>%
    t() %>%
    write.table(file = "processed/STEM_samples_info_old_sigma.tsv", sep = "\t",
                row.names = TRUE, quote = FALSE, col.names = FALSE, na = "")
meta_o %>%
    dplyr::filter(colname %in% colnames(diff)[-1]) %>%
    dplyr::select(-colname) %>%
    t() %>%
    write.table(file = "processed/DIFF_samples_info_old_sigma.tsv", sep = "\t",
                row.names = TRUE, quote = FALSE, col.names = FALSE, na = "")

val <- sapply(meta_o[2:ncol(meta_o)], unique)

data.frame(pos = "SAMPLEINFO",
           value = rep(names(val), lengths(val)),
           "s" = unlist(val, use.names = FALSE),
           color = "blue") %>%
    mutate(color = is.character(color)) %>%
    mutate(color = case_when(s == "inflamed" ~ "red",
                             s == "healthy" ~ "green",
                             s == "not involved" ~ "LIGHT_YELLOW",
                             s == "not inflamed" ~ "orange",
                             s == "CTRL" ~ "yellow",
                             s == "CD" ~ "blue",
                             s == "pediatric" ~ "CYAN",
                             s == "adult" ~ "light_orange",
                             s == "sigma" ~ "white",
                             s == "ileum" ~ "light_blue"
    )) %>%
    filter(!is.na(color)) %>%
    mutate(color = toupper(color)) %>%
    write.table(file = "processed/colors_info_old_sigma.tsv", quote = FALSE, sep = "\t",
                col.names = FALSE, row.names = FALSE)

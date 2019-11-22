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


# Design ####
grouping <- apply(meta[, c("cell_type", "GROUP", "TYPE", "LOCATION"), drop = FALSE], 1,
                  paste0, collapse = "_")
grouping <- as.factor(grouping)
design <- model.matrix(~0 + grouping, grouping)
colnames(design) <- levels(grouping)
design <- cbind(Intercept = 1, design) # Truco 1 afegir intercept
rownames(design) <- meta$colname

design <- design[ , colSums(design) != 1] # Truco 2 treure aquelles amb 1 sola mostra

# Add id of the sample
# samples_d <- model.matrix(~0 + SAMPLE2, data = meta[, "SAMPLE2", drop = FALSE])
# design <- cbind(design, samples_d)
# design <- design[ , colSums(design) != 1]

# Contrasts ####
# Helper functions for contrasts
paste_plus <- function(...){paste(..., collapse = " + ")}
paste_ab <- function(a, b){paste(a, "- (", b, ")")}
grep_level <- function(x, a){
    # Force to actually look up on the design being used!
    # This makes it dependent on the order it is run!!
    grep(x, a, value = TRUE)
}
contrasting <- function(x, y, design_m) {
    d <- colnames(design_m)
    paste_ab(paste_plus(grep_level(x, d)),
             paste_plus(grep_level(y, d)))
}

# DIFF eq 1, STEM to -1
contrasts <- c(
    diff_stem = contrasting("^DIFF", "^STEM", design), # Siempre por separado
    CD_CTRL = contrasting("_CD_", "_CTRL_", design),
    STEM__ileum_sigma = contrasting("STEM.*ileum", "STEM.*sigma", design), # STEM ileum vs STEM sigma
    DIFF__ileum_sigma = contrasting("DIFF.*ileum", "DIFF.*sigma", design),# DIFF ileum vs DIFF sigma
    pediatric_adult = contrasting("pediatric", "adult", design),
    STEM__pediatric_adult = contrasting("STEM.*pediatric", "STEM.*adult", design),
    DIFF__pediatric_adult = contrasting("DIFF.*pediatric", "DIFF.*adult", design),
    ileum_sigma = contrasting("ileum", "sigma", design),
    ileum__STEM_DIFF = contrasting("STEM.*ileum", "DIFF.*ileum", design), # OK
    sigma__STEM_DIFF = contrasting("STEM.*sigma", "DIFF.*sigma", design), # OK
    sigma_DIFF__pediatric_adult = contrasting("DIFF.*pediatric_sigma", "DIFF.*adult_sigma", design),
    sigma_STEM__pediatric_adult = contrasting("STEM.*pediatric_sigma", "STEM.*adult_sigma", design),
    sigma__pediatric_adult = contrasting("pediatric_sigma", "adult_sigma", design),
    ileum__pediatric_adult = contrasting("pediatric_ileum", "adult_ileum", design),
    sigma_diff__pediatric_adult = contrasting("DIFF.*pediatric_sigma", "DIFF.*adult_sigma", design),
    STEM_ileum_CD__pediatric_adult = contrasting("STEM_CD_pediatric_ileum", "STEM_CD_adult_ileum", design),
    STEM_ileum__pediatric_adult = contrasting("STEM_.*_pediatric_ileum", "STEM_.*_adult_ileum", design),
    STEM_sigma__pediatric_adult = contrasting("STEM_.*_pediatric_sigma", "STEM_.*_adult_sigma", design),
    DIFF_sigma__pediatric_adult = contrasting("DIFF_.*_pediatric_sigma", "DIFF_.*_adult_sigma", design)
)

contr.matrix <- makeContrasts(contrasts = contrasts, levels = colnames(design))
colnames(contr.matrix) <- names(contrasts)


pick <- design %*% contr.matrix # To know which samples are taken
# Look if there is any contrast with just one sample on one side
p <- t(apply(pick, 2, table))

dge_rna <- DGEList(expr)
keep <- filterByExpr(dge_rna)
dge_rna <- dge_rna[keep, ]

# Prepare files for Juanjo
write.csv(dge_rna$counts, file = "processed/old_counts.csv",
          row.names = TRUE)
d2 <- pick
d2[d2 == -1] <- 2
d2[d2 == 0] <- NA
write.csv(d2, file = "processed/comparisons.csv", row.names = TRUE)

dge_rna <- calcNormFactors(dge_rna, method = "TMM")
# Limma p 125-126 two vooms
voom_rna <- voom(dge_rna, plot = TRUE, normalize.method = "cyclicloess")

vfit <- lmFit(voom_rna, design = design)

# Continue evaluation whole ####
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
ttc <- vector("list", 21)
pdf("processed/MA_complex_old_samples2.pdf")
for (i in seq_len(ncol(contr.matrix))) {
    ttc[[i]] <- topTable(efit, coef = names(contrasts)[i],
                         number = Inf)
    plot(ttc[[i]]$AveExpr, ttc[[i]]$logFC,
         main = names(contrasts)[i], pch = ".")
}
dev.off()


dt <- decideTests(efit, lfc = log2(1.5), adjust.method = "fdr")
dtt <- summary(dt)
dtt

# Evaluation ####



tt <- vector("list", length(contrasts))
pdf("processed/MA_single.pdf")
for (i in seq_along(contrasts)) {
# Checking contrasts  (order of contrast affects them )
# they are picking the right samples

    contr <- pick[, i]
    dge_rna <- DGEList(expr[, contr != 0])
    keep <- filterByExpr(dge_rna)
    dge_rna <- dge_rna[keep, , keep.lib.size = FALSE]
    dge_rna <- calcNormFactors(dge_rna, method = "TMM")
    logCPM <- cpm(dge_rna, log=TRUE, prior.count=3)
    voom_rna <- voom(logCPM, plot = FALSE, normalize.method = "quantile")
    stopifnot(ncol(voom_rna$E) == p[i])
    boxplot(voom_rna$E, main = names(contrasts)[i],
            col = as.factor(contr[contr != 0]))

    vfit <- lmFit(voom_rna, contr[contr != 0])
    efit <- eBayes(vfit)

    tt[[i]] <- topTreat(efit, coef = 1, number = Inf)
    plot(tt[[i]]$AveExpr, tt[[i]]$logFC, main = names(contrasts)[i], pch = 16)

}
dev.off()

# Prepare for GET ####
# perl ~/Documents/projects/GETS/gets.pl --matrix=./processed/STEM_old.tsv --geneinfo=./processed/gene_info_old.tsv --sampleinfo=./processed/STEM_samples_info_old.tsv --colors=./processed/colors_info_old.tsv --output=./processed/STEM_GETS_old --center=TRUE --overwrite=TRUE
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


write.table(gene_info2[, -1], "processed/gene_info_old.tsv", sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE, na = "")

stopifnot(ncol(voom_rna$E) == nrow(meta))
r <- order(meta$TYPE, meta$cell_type, meta$LOCATION, meta$GROUP)
meta_o <- meta[r, -10]
meta_o <- meta_o[, c("colname", "TYPE", "cell_type", "LOCATION", "GROUP", "status",
           "age", "AaD", "SAMPLE", "SAMPLE2")]
expr_csv <- voom_rna$E[, r]
stopifnot(ncol(expr_csv) == nrow(meta))
stopifnot(nrow(gene_info2) == nrow(expr_csv))
expr2_csv <- cbind.data.frame(genes = rownames(expr_csv), expr_csv)
write.table(expr2_csv, file = "processed/old.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
stem <- cbind.data.frame(genes = rownames(expr_csv), expr_csv[, meta_o$cell_type == "STEM"])
write.table(stem,
            file = "processed/STEM_old.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")
diff <- cbind.data.frame(genes = rownames(expr_csv), expr_csv[, meta_o$cell_type == "DIFF"])
write.table(diff,
            file = "processed/DIFF_old.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")

meta_o %>%
    dplyr::filter(colname %in% colnames(expr2_csv)) %>%
    dplyr::select(-colname) %>%
    t() %>%
    write.table(file = "processed/samples_info_old.tsv", sep = "\t",
                row.names = TRUE, quote = FALSE, col.names = FALSE, na = "")
meta_o %>%
    dplyr::filter(colname %in% colnames(stem)[-1]) %>%
    dplyr::select(-colname) %>%
    t() %>%
    write.table(file = "processed/STEM_samples_info_old.tsv", sep = "\t",
                row.names = TRUE, quote = FALSE, col.names = FALSE, na = "")
meta_o %>%
    dplyr::filter(colname %in% colnames(diff)[-1]) %>%
    dplyr::select(-colname) %>%
    t() %>%
    write.table(file = "processed/DIFF_samples_info_old.tsv", sep = "\t",
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
    write.table(file = "processed/colors_info_old.tsv", quote = FALSE, sep = "\t",
                col.names = FALSE, row.names = FALSE)

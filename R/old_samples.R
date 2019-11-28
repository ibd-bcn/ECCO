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
shared_samples <- intersect(colnames(umat), meta$colname[meta$reanalyzed])
expr <- umat[, shared_samples]
meta <- meta[match(shared_samples, meta$colname), ]

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
# On 25/11/2019 Isa asked me to add all the comparisons instead of just the few
# She also identified some duplicated comparisons (my fault)
# On 28/11/2019 morning Isa asked me to make it as in the previous analysis of microarrays
# On 28/112019 afternon Isa asked me to make it the full 36 contrasts from pegs_organoids_2017 project
# and remove those that don't make sense
contrasts <- c(
    stem_diff_colon_A = contrasting("STEM.*adult_sigma", "DIFF.*adult_sigma", design), # 1
    stem_diff_colon_P = contrasting("STEM.*pediatric_sigma", "DIFF.*pediatric_sigma", design), #  2
    stem_diff_ileum_A = contrasting("STEM.*adult_ileum", "DIFF.*adult_ileum", design), #  3
    # stem_diff_ileum_P = contrasting("STEM.*pediatric_ileum", "DIFF.*pediatric_ileum", design),  #  4
    stem_ileum_colon_A = contrasting("STEM.*adult_ileum", "STEM.*adult_sigma", design), #  5
    stem_ileum_colon_P = contrasting("STEM.*pediatric_ileum", "STEM.*pediatric_sigma", design), #  6
    diff_ileum_colon_A = contrasting("DIFF.*adult_ileum", "DIFF.*adult_sigma", design), #  7
    # diff_ileum_colon_P = contrasting("DIFF.*pediatric_ileum", "DIFF.*pediatric_sigma", design),  #  8
    diff_H_CD_colon_A = contrasting("DIFF_CTRL_adult_sigma", "DIFF_CD_adult_sigma", design), #  9
    stem_H_CD_colon_A = contrasting("STEM_CTRL_adult_sigma", "STEM_CD_adult_sigma", design),  #  10
    stem_H_CD_colon_P = contrasting("STEM_CTRL_pediatric_sigma", "STEM_CD_pediatric_sigma", design), #  11
    # diff_H_CD_colon_P = contrasting("DIFF_CTRL_pediatric_sigma", "DIFF_CD_pediatric_sigma", design), #  12
    stem_ileum_CD_A_P = contrasting("STEM_CD_adult_ileum", "STEM_CD_pediatric_ileum", design), #  13
    stem_colon_CD_A_P = contrasting("STEM_CD_adult_sigma", "STEM_CD_pediatric_sigma", design), #  14
    # diff_ileum_CD_A_P = contrasting("DIFF_CD_adult_ileum", "DIFF_CD_pediatric_ileum", design), #  15
    diff_colon_CD_A_P = contrasting("DIFF_CD_adult_sigma", "DIFF_CD_pediatric_sigma", design), #  16
    #  17
    #  18
    #  19
    #  20
    #  21
    #  22
    #  23
    #  24
    #  25
    #  26
    #  27
    diff_colon_H_CD = contrasting("DIFF_CTRL_.*_sigma", "DIFF_CD_.*_sigma", design), #  28
    stem_colon_CD_A_P = contrasting("STEM_CD_adult_sigma", "STEM_CD_pediatric_sigma", design), #  29
    diff_colon_CD_A_P = contrasting("DIFF_CD_adult_sigma", "DIFF_CD_pediatric_sigma", design), #  30
    stem_ileum_CD_A_P = contrasting("STEM_CD_adult_ileum", "STEM_CD_pediatric_ileum", design), #  31
    # DIFF_ileum_CD_A_P = contrasting("DIFF_CD_adult_ileum", "DIFF_CD_pediatric_ileum", design), #  32
    stem_ileum_P_A = contrasting("STEM.*pediatric_ileum", "STEM.*adult_ileum", design), #  33
    stem_colon_P_A = contrasting("STEM.*pediatric_sigma", "STEM.*adult_sigma", design), #  34
    # diff_ileum_P_A = contrasting("DIFF.*pediatric_ileum", "DIFF*adult_ileum", design), #  35 # Impossible
    diff_colon_P_A = contrasting("DIFF.*pediatric_sigma", "DIFF.*adult_sigma", design), #  36
    diff_vs_stem = contrasting("^DIFF", "^STEM", design), # Siempre por separado
    CD_vs_CTRL = contrasting("_CD_", "_CTRL_", design),
    pediatric_vs_adult = contrasting("pediatric", "adult", design),
    ileum_vs_sigma = contrasting("ileum", "sigma", design)

    # STEM__ileum_vs_sigma = contrasting("STEM.*ileum", "STEM.*sigma", design), # STEM ileum vs STEM sigma
    # DIFF__ileum_vs_sigma = contrasting("DIFF.*ileum", "DIFF.*sigma", design),# DIFF ileum vs DIFF sigma
    # STEM__pediatric_vs_adult = contrasting("STEM.*pediatric", "STEM.*adult", design),
    # DIFF__pediatric_vs_adult = contrasting("DIFF.*pediatric", "DIFF.*adult", design),
    # ileum__STEM_vs_DIFF = contrasting("STEM.*ileum", "DIFF.*ileum", design), # OK
    # sigma__STEM_vs_DIFF = contrasting("STEM.*sigma", "DIFF.*sigma", design), # OK
    # # sigma_DIFF__pediatric_vs_adult = contrasting("DIFF.*pediatric_sigma", "DIFF.*adult_sigma", design),
    # # sigma_STEM__pediatric_vs_adult = contrasting("STEM.*pediatric_sigma", "STEM.*adult_sigma", design),
    # sigma__pediatric_vs_adult = contrasting("pediatric_sigma", "adult_sigma", design),
    # ileum__pediatric_vs_adult = contrasting("pediatric_ileum", "adult_ileum", design),
    # # STEM_ileum_CD__pediatric_vs_adult = contrasting("STEM_CD_pediatric_ileum", "STEM_CD_adult_ileum", design),
    # # STEM_ileum__pediatric_adult = contrasting("STEM_.*_pediatric_ileum", "STEM_.*_adult_ileum", design),
    # # STEM_CD_ileum__pediatric_adult = contrasting("STEM_CD_pediatric_ileum", "STEM_CD_adult_ileum", design),
    # # DIFF_CD_ileum__pediatric_adult = contrasting("DIFF_CD_pediatric_ileum", "DIFF_CD_adult_ileum", design),
    # DIFF_CD_sigma__pediatric_vs_adult = contrasting("DIFF_CD_pediatric_sigma", "DIFF_CD_adult_sigma", design),
    # STEM_CD_sigma__pediatric_vs_adult = contrasting("STEM_CD_pediatric_sigma", "STEM_CD_adult_sigma", design),
    #
    # # adult_CD_ileum__STEM_vs_DIFF = contrasting("STEM_CD_adult_ileum", "DIFF_CD_adult_ileum", design),
    # # pediatric_CD_ileum__STEM_DIFF = contrasting("STEM_CD_pediatric_ileum", "DIFF_CD_pediatric_ileum", design),
    # # CD_ileum__STEM_DIFF = contrasting("STEM_CD_.*_ileum", "DIFF_CD_.*_ileum", design),
    # CD_sigma__STEM_vs_DIFF = contrasting("STEM_CD_.*_sigma", "DIFF_CD_.*_sigma", design),
    #
    # # adult_CTRL_ileum__STEM_DIFF = contrasting("STEM_CTRL_adult_ileum", "DIFF_CTRL_adult_ileum", design),
    # # pediatric_CTRL_ileum__STEM_DIFF = contrasting("STEM_CTRL_pediatric_ileum", "DIFF_CTRL_pediatric_ileum", design),
    # # CTRL_ileum__STEM_DIFF = contrasting("STEM_CTRL_.*_ileum", "DIFF_CTRL_.*_ileum", design),
    # CTRL_sigma__STEM_vs_DIFF = contrasting("STEM_CTRL_.*_sigma", "DIFF_CTRL_.*_sigma", design),
    # # CD_sigma_DIFF__pediatric_adult = contrasting("DIFF_CD_pediatric_sigma", "DIFF_CD_adult_sigma", design),
    # CTRL_sigma_STEM__pediatric_vs_adult = contrasting("STEM_CTRL_pediatric_sigma", "STEM_CTRL_adult_sigma", design),
    #
    # STEM_sigma_pediatric__CD_vs_CTRL = contrasting("STEM_CD_pediatric_sigma", "STEM_CTRL_pediatric_sigma", design),
    # STEM_sigma_adult__CD_vs_CTRL = contrasting("STEM_CD_adult_sigma", "STEM_CTRL_adult_sigma", design),
    # # DIFF_sigma_pediatric_CD_CTRL = contrasting("DIFF_CD_pediatric_sigma", "DIFF_CTRL_pediatric_sigma", design),
    # DIFF_sigma_adult__CD_vs_CTRL = contrasting("DIFF_CD_adult_sigma", "DIFF_CTRL_adult_sigma", design),
    # DIFF_CD_adult__ileum_vs_sigma = contrasting("DIFF_CD_adult_ileum", "DIFF_CD_adult_sigma", design),
    # DIFF_CD__adult_ileum_vs_pediatric_sigma = contrasting("DIFF_CD_adult_ileum", "DIFF_CD_pediatric_sigma", design),
    # DIFF_adult__CD_ileum_vs_CTRL_sigma = contrasting("DIFF_CD_adult_ileum", "DIFF_CTRL_adult_sigma", design),
    # # CD_adult_ileum__DIFF_vs_STEM = contrasting("DIFF_CD_adult_ileum", "STEM_CD_adult_ileum", design),
    # CD_adult__DIFF_ileum_vs_STEM_sigma = contrasting("DIFF_CD_adult_ileum", "STEM_CD_adult_sigma", design),
    # CD__DIFF_adult_ileum_vs_STEM_pediatric_ileum = contrasting("DIFF_CD_adult_ileum", "STEM_CD_pediatric_ileum", design),
    # CD__DIFF_adult_ileum_vs_STEM_pediatric_sigma = contrasting("DIFF_CD_adult_ileum", "STEM_CD_pediatric_sigma", design),
    # adult__DIFF_CD_ileum_vs_STEM_CTRL_sigma = contrasting("DIFF_CD_adult_ileum", "STEM_CTRL_adult_sigma", design),
    # DIFF_CD_adult_ileum_vs_STEM_CTRL_pediatric_sigma = contrasting("DIFF_CD_adult_ileum", "STEM_CTRL_pediatric_sigma", design),
    # DIFF_adult_sigma__CD_vs_CTRL = contrasting("DIFF_CD_adult_sigma", "STEM_CD_adult_ileum", design),
    # CD_adult_sigma__DIFF_vs_STEM = contrasting("DIFF_CD_adult_sigma", "STEM_CD_adult_sigma", design),
    # CD__DIFF_adult_sigma_vs_STEM_pediatric_ileum = contrasting("DIFF_CD_adult_sigma", "STEM_CD_pediatric_ileum", design),
    # CD_sigma__DIFF_adult_vs_STEM_pediatric = contrasting("DIFF_CD_adult_sigma", "STEM_CD_pediatric_sigma", design),
    # adult_sigma__DIFF_CD_vs_STEM_CTRL = contrasting("DIFF_CD_adult_sigma", "STEM_CTRL_adult_sigma", design),
    # sigma__DIFF_CD_adult_vs_STEM_CTRL_pediatric = contrasting("DIFF_CD_adult_sigma", "STEM_CTRL_pediatric_sigma", design),
    # DIFF_sigma__CD_pediatric_vs_CTRL_adult = contrasting("DIFF_CD_pediatric_sigma", "DIFF_CTRL_adult_sigma", design),
    # CD__DIFF_pediatric_sigma_vs_STEM_adult_ileum = contrasting("DIFF_CD_pediatric_sigma", "STEM_CD_adult_ileum", design),
    # CD_sigma__DIFF_CD_vs_STEM_CTRL = contrasting("DIFF_CD_pediatric_sigma", "STEM_CTRL_pediatric_sigma", design),
    # # CD_pediatric__DIFF_sigma_STEM_ileum = contrasting("DIFF_CD_adult_sigma", "STEM_CD_pediatric_ileum", design),
    # # 22 of combn(colnames(design)[-1], 2):
    # CD_pediatric_sigma__DIFF_vs_STEM = contrasting("DIFF_CD_pediatric_sigma", "STEM_CD_pediatric_sigma", design),
    # sigma__DIFF_CD_pediatric_vs_STEM_CTRL_adult_sigma = contrasting("DIFF_CD_pediatric_sigma", "STEM_CTRL_adult_sigma", design),
    # # pediatric_sigma__DIFF_CD_STEM_CTRL = contrasting("DIFF_CD_pediatric_sigma", "STEM_CTRL_pediatric_sigma", design),
    # adult__DIFF_CTRL_sigma_vs_STEM_CD_ileum = contrasting("DIFF_CTRL_adult_sigma", "STEM_CD_adult_ileum", design),
    # adult_sigma__DIFF_CTRL_vs_STEM_CD = contrasting("DIFF_CTRL_adult_sigma", "STEM_CD_adult_sigma", design),
    # # 27
    # DIFF_CTRL_adult_sigma_vs_STEM_CD_adult_sigma = contrasting("DIFF_CTRL_adult_sigma", "STEM_CD_pediatric_ileum", design),
    # sigma__DIFF_CTRL_adult_vs_STEM_CD_pediatric = contrasting("DIFF_CTRL_adult_sigma", "STEM_CD_pediatric_sigma", design),
    # CTRL_adult_sigma__DIFF_vs_STEM = contrasting("DIFF_CTRL_adult_sigma", "STEM_CTRL_adult_sigma", design),
    # CTRL_sigma__DIFF_adult_vs_STEM_pediatric = contrasting("DIFF_CTRL_adult_sigma", "STEM_CTRL_pediatric_sigma", design),
    # STEM_CD_adult__ileum_vs_sigma = contrasting("STEM_CD_adult_ileum", "STEM_CD_adult_sigma", design),
    # # 32
    # # STEM_CD_ileum__adult_vs_pediatric = contrasting("STEM_CD_adult_ileum", "STEM_CD_pediatric_ileum", design),
    # STEM_CD__adult_ileum_vs_pediatric_sigma = contrasting("STEM_CD_adult_ileum", "STEM_CD_pediatric_sigma", design),
    # STEM_adult__CD_ileum_vs_CTRL_sigma = contrasting("STEM_CD_adult_ileum", "STEM_CTRL_adult_sigma", design),
    # STE__CD_adult_ileum_vs_CTRL_pediatric_sigma = contrasting("STEM_CD_adult_ileum", "STEM_CTRL_pediatric_sigma", design),
    # STEM_CD__adult_sigma_vs_pediatric_ileum = contrasting("STEM_CD_adult_sigma", "STEM_CD_pediatric_ileum", design),
    # # 37
    # # STEM_CD_sigma__adult_pediatric = contrasting()
    # # STEM_sigma_adult__CD_CTRL
    # STEM_sigma__CD_adult_vs_CTRL_pediatric = contrasting("STEM_CD_adult_sigma", "STEM_CTRL_pediatric_sigma", design),
    # STEM_CD_pediatric__ileum_vs_sigma = contrasting("STEM_CD_pediatric_ileum", "STEM_CD_pediatric_sigma", design),
    # STEM__CD_pediatric_ileum_vs_CTRL_adult_sigma = contrasting("STEM_CD_pediatric_ileum", "STEM_CTRL_adult_sigma", design),
    # # 42
    # STEM_pediatric__CD_ileum_vs_CTRL_sigma = contrasting("STEM_CD_pediatric_ileum", "STEM_CTRL_pediatric_sigma", design),
    # STEM_sigma__CD_pediatric_vs_CTRL_adult = contrasting("STEM_CD_pediatric_sigma", "STEM_CTRL_adult_sigma", design),
    # # STEM_pediatric_sigma_CD_CTRL = contrasting("STEM_CD_pediatric_sigma", "STEM_CTRL_pediatric_sigma", design),
    # STEM_CTRL_sigma__adult_vs_pediatric = contrasting("STEM_CTRL_adult_sigma", "STEM_CTRL_pediatric_sigma", design)
    )

contr.matrix <- makeContrasts(contrasts = contrasts, levels = colnames(design))

design2 <- meta %>%
    mutate(diff_colon_CD_Healthy_involved_A = case_when(
        cell_type == "DIFF" & status == "healthy" & GROUP == "CD" & TYPE == "adult" & LOCATION == "sigma"~ 1,
        cell_type == "DIFF" & status == "involved" & GROUP == "CD" & TYPE == "adult" & LOCATION == "sigma" ~ -1,
        TRUE ~ 0),
        ) %>%
    select(diff_colon_CD_Healthy_involved_A,
           ) %>%
    as.matrix()
contr.matrix <- makeContrast(contrasts = contrasts2, levels = colnames(design2))
colnames(contr.matrix) <- names(contrasts)

pick <- design %*% contr.matrix # To know which samples are taken

samples_used <- t(apply(pick, 2, table))
stopifnot(ncol(samples_used) == 3)

pu <- unique(t(pick))
# If there are duplicate contrast check the output and remove them
if (nrow(t(pick)) - nrow(pu) != 0){

    duplicated_contr <- setdiff(rownames(t(pick)), rownames(pu))

    dups <- sapply(duplicated_contr, function(y) {
        apply(pick, 2, function(x) {all(x ==pick[, y])})
    })
    dups
}

# Look if there is any contrast with just one sample on one side
dge_rna <- DGEList(expr)
keep <- filterByExpr(dge_rna)
dge_rna <- dge_rna[keep, ]
CPM <- cpm(dge_rna$counts, log = TRUE)
keep <- rowSums(CPM > 0.5) >= 23 # 24 samples with more expression than 1 CPM
keep2 <- rowSums(CPM)
# Prepare files for Juanjo
write.csv(dge_rna$counts, file = "processed/old_counts.csv",
          row.names = TRUE)

d2 <- pick
d2[d2 == 1] <- 2
d2[d2 == -1] <- 1
d2[d2 == 0] <- NA
write.csv(d2, file = "processed/comparisons.csv", row.names = TRUE)

dge_rna <- calcNormFactors(dge_rna, method = "TMM")
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

saveRDS(ttc, "processed/tt_complex.RDS")

dt <- decideTests(efit, lfc = log2(1.5), adjust.method = "fdr")
dtt <- summary(dt)
dtt

# Evaluation ####
tt <- vector("list", length(contrasts))

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

r <- order(meta$cell_type, meta$LOCATION, meta$TYPE, meta$GROUP)
meta_o <- meta[r, -10]
meta_o <- meta_o[, c("colname", "cell_type", "LOCATION", "TYPE", "GROUP", "status",
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

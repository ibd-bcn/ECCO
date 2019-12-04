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

# Subset to just new samples
expr <- umat[, colnames(umat) %in% meta$colname[!meta$reanalyzed]]
meta <- meta[meta$colname %in% colnames(expr), ]

dge_rna <- DGEList(expr)
keep <- filterByExpr(dge_rna)
dge_rna <- dge_rna[keep, ]
dge_rna <- calcNormFactors(dge_rna)
voom_rna <- voom(dge_rna, plot = FALSE, normalize.method = "quantile")

# Design ####
grouping <- apply(meta[, c("cell_type", "GROUP", "TYPE", "LOCATION")], 1,
                  paste0, collapse = "_")
grouping <- as.factor(grouping)
design <- model.matrix(~0 + grouping, grouping)
colnames(design) <- levels(grouping)
# design <- cbind(design, Reanalyzed = as.numeric(meta$reanalyzed))

# corfit <- duplicateCorrelation(dge_rna$counts, design, block = meta$SAMPLE2)
# samples_d <- model.matrix(~0 + SAMPLE2, data = meta[, "SAMPLE2", drop = FALSE])
# design <- cbind(design, samples_d)

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

# On 25/11/2019 Isa asked me to add all the comparisons instead of just the few
# She also identified some duplicated comparisons (my fault)
# On 28/11/2019 morning Isa asked me to make it as in the previous analysis of microarrays
# On 28/112019 afternon Isa asked me to make it the full 36 contrasts from pegs_organoids_2017 project
# and remove those that don't make sense
# On 29/11/2019 afternon I asked Isa  about what to do about contrast14_stem_colon_CD_A_P and contrast29_CD_A_P
# The difference is between having healthy and not having healthy samples on each group
contrasts <- c(
    c1_stem_diff_colon_A = contrasting("STEM.*adult_sigma", "DIFF.*adult_sigma", design), # 1
    c2_stem_diff_colon_P = contrasting("STEM.*pediatric_sigma", "DIFF.*pediatric_sigma", design), #  2
    c3_stem_diff_ileum_A = contrasting("STEM.*adult_ileum", "DIFF.*adult_ileum", design), #  3
    c4_stem_diff_ileum_P = contrasting("STEM.*pediatric_ileum", "DIFF.*pediatric_ileum", design),  #  4
    c5_stem_ileum_colon_A = contrasting("STEM.*adult_ileum", "STEM.*adult_sigma", design), #  5
    c6_stem_ileum_colon_P = contrasting("STEM.*pediatric_ileum", "STEM.*pediatric_sigma", design), #  6
    c7_diff_ileum_colon_A = contrasting("DIFF.*adult_ileum", "DIFF.*adult_sigma", design), #  7
    c8_diff_ileum_colon_P = contrasting("DIFF.*pediatric_ileum", "DIFF.*pediatric_ileum", design),  #  8
    # c9_diff_H_CD_colon_A in reality has not inflamed CD only!! vs healthy ctrl
    # c10_stem_H_CD_colon_A in reality has non inflamed CD only!! vs healthy ctrl
    c11_stem_H_CD_colon_P = contrasting("STEM_CTRL_pediatric_sigma", "STEM_CD_pediatric_sigma", design), #  11
    c12_diff_H_CD_colon_P = contrasting("DIFF_CTRL_pediatric_sigma", "DIFF_CD_pediatric_sigma", design), #  12
    # c13_stem_ileum_CD_A_P in reality stem_ileum_CD_A_P has not inflamed only!!
    # c14_stem_colon_CD_A_P has not inflamed only on adults!!
    # c15_diff_ileum_CD_A_P only not inflamed!!
    # c16_diff_colon_CD_A_P adult diff not inflamed!!
    # 17
    # 18
    # 19
    # 20
    # 21
    # 22
    # 23
    # 24
    # 25
    # 26
    # 27
    # 28
    c29_stem_colon_CD_A_P = contrasting("STEM_CD_adult_sigma", "STEM_CD_pediatric_sigma", design), #  29
    c30_diff_colon_CD_A_P = contrasting("DIFF_CD_adult_sigma", "DIFF_CD_pediatric_sigma", design), #  30
    c31_stem_ileum_CD_A_P = contrasting("STEM_CD_adult_ileum", "STEM_CD_pediatric_ileum", design), #  31
    c32_diff_ileum_CD_A_P = contrasting("DIFF_CD_adult_ileum", "DIFF_CD_pediatric_ileum", design), # 32 No n
    # Just in case Siempre por separado
    c33_diff_vs_stem = contrasting("^DIFF", "^STEM", design),
    c34_CD_vs_CTRL = contrasting("_CD_", "_CTRL_", design),
    c35_pediatric_vs_adult = contrasting("pediatric", "adult", design),
    c36_ileum_vs_sigma = contrasting("ileum", "sigma", design)
)

contr.matrix <- makeContrasts(contrasts = contrasts, levels = colnames(design))
colnames(contr.matrix) <- names(contrasts)

design2 <- meta %>%
    mutate(
        c9_diff_H_CD_colon_A = case_when( # Number 9
            cell_type == "DIFF" & status == "not inflamed" & GROUP == "CD" & TYPE == "adult" & LOCATION == "sigma" ~ -1,
            cell_type == "DIFF" & status == "healthy" & GROUP == "CTRL" & TYPE == "adult" & LOCATION == "sigma"~ 1,
            TRUE ~ 0),
        c10_stem_H_CD_colon_A = case_when( # Number 10
            cell_type == "STEM" & status == "not inflamed" & GROUP == "CD" & TYPE == "adult" & LOCATION == "sigma" ~ -1,
            cell_type == "STEM" & status == "healthy" & GROUP == "CTRL" & TYPE == "adult" & LOCATION == "sigma"~ 1,
            TRUE ~ 0),
        c13_stem_ileum_CD_A_P = case_when( # Number 13
            cell_type == "STEM" & status == "not inflamed" & GROUP == "CD" & TYPE == "pediatric" & LOCATION == "ileum" ~ -1,
            cell_type == "STEM" & status == "not inflamed" & GROUP == "CD" & TYPE == "adult" & LOCATION == "ileum"~ 1,
            TRUE ~ 0),
        c14_stem_colon_CD_A_P = case_when( # Number 14
            cell_type == "STEM" & GROUP == "CD" & TYPE == "pediatric" & LOCATION == "sigma" ~ -1,
            cell_type == "STEM" & status == "not inflamed" & GROUP == "CD" & TYPE == "pediatric" & LOCATION == "sigma"~ 1,
            TRUE ~ 0),
        c15_diff_ileum_CD_A_P = case_when( # Number 15
            cell_type == "DIFF" & status == "not inflamed" & GROUP == "CD" & TYPE == "pediatric" & LOCATION == "ileum" ~ -1,
            cell_type == "DIFF" & status == "not inflamed" & GROUP == "CD" & TYPE == "adult" & LOCATION == "ileum"~ 1,
            TRUE ~ 0),
        c16_diff_colon_CD_A_P = case_when( # Number 16
            cell_type == "DIFF" & GROUP == "CD" & TYPE == "pediatric" & LOCATION == "sigma" ~ -1,
            cell_type == "DIFF" & status == "not inflamed" & GROUP == "CD" & TYPE == "adult" & LOCATION == "sigma"~ 1,
            TRUE ~ 0),
        c17_colon_CD_healthy_involved_A = case_when( # Number 17
            cell_type == "DIFF" & status == "not inflamed" & GROUP == "CD" & TYPE == "adult" & LOCATION == "sigma" ~ -1,
            cell_type == "DIFF" & status == "healthy" & GROUP == "CD" & TYPE == "adult" & LOCATION == "sigma"~ 1,
            TRUE ~ 0),
        c18_stem_colon_CD_Healthy_involved_A = case_when( # Number 18
            cell_type == "STEM" & status == "healthy" & GROUP == "CTRL" & TYPE == "adult" & LOCATION == "sigma"~ 1,
            cell_type == "STEM" & status == "not inflamed" & GROUP == "CD" & TYPE == "adult" & LOCATION == "sigma" ~ -1,
            TRUE ~ 0),
        c19_diff_ileum_CD_healthy_involved_A = case_when( # Number 19
            cell_type == "DIFF" & status == "healthy" & GROUP == "CTRL" & TYPE == "adult" & LOCATION == "ileum"~ 1,
            cell_type == "DIFF" & status == "not inflamed" & GROUP == "CD" & TYPE == "adult" & LOCATION == "ileum" ~ -1,
            TRUE ~ 0),
        c20_stem_ileum_CD_healthy_involved_A = case_when( # Number 20
            cell_type == "STEM" & status == "healthy" & GROUP == "CTRL" & TYPE == "adult" & LOCATION == "ileum"~ 1,
            cell_type == "STEM" & status == "not inflamed" & GROUP == "CD" & TYPE == "adult" & LOCATION == "ileum" ~ -1,
            TRUE ~ 0),
        c21_diff_colon_CD_healthy_involved_P = case_when( # Number 21
            cell_type == "DIFF" & status == "healthy" & GROUP == "CC" & TYPE == "pediatric" & LOCATION == "sigma"~ 1,
            cell_type == "DIFF" & GROUP == "CD" & TYPE == "pediatric" & LOCATION == "sigma" ~ -1,
            TRUE ~ 0),
        c22_stem_colon_CD_Healthy_involved_P = case_when( # Number 22
            cell_type == "STEM" & status == "healthy" & GROUP == "CTRL" & TYPE == "pediatric" & LOCATION == "sigma"~ 1,
            cell_type == "STEM" & status == "not inflamed" & GROUP == "CD" & TYPE == "pediatric" & LOCATION == "sigma" ~ -1,
            TRUE ~ 0),
        c23_diff_ileum_CD_healthy_involved_P = case_when( # Number 23
            cell_type == "DIFF" & status == "healthy" & GROUP == "CTRL" & TYPE == "pediatric" & LOCATION == "ileum"~ 1,
            cell_type == "DIFF" & status == "not inflamed" & GROUP == "CD" & TYPE == "pediatric" & LOCATION == "ileum" ~ -1,
            TRUE ~ 0),
        c24_stem_ileum_CD_healthy_involved_P = case_when( # Number 24
            cell_type == "STEM" & status == "healthy" & GROUP == "CTRL" & TYPE == "pediatric" & LOCATION == "ileum"~ 1,
            cell_type == "STEM" & status == "involved" & GROUP == "CD" & TYPE == "pediatric" & LOCATION == "ileum" ~ -1,
            TRUE ~ 0),
        c25_stem_ileum_H_healthy_CD = case_when( # Number 25 just 1 sample
            cell_type == "STEM" & GROUP == "CTRL" & LOCATION == "ileum"~ 1,
            cell_type == "STEM" & status == "not involved" & GROUP == "CD" & LOCATION == "ileum" ~ -1,
            TRUE ~ 0),
        c26_diff_ileum_H_healthy_CD = case_when( # Number 26
            cell_type == "DIFF" & GROUP == "CTRL" & LOCATION == "ileum"~ 1,
            cell_type == "DIFF" & status == "not involved" & GROUP == "CD" & LOCATION == "ileum" ~ -1,
            TRUE ~ 0),
        c27_stem_colon_H_healthy_CD = case_when( # Number 27
            cell_type == "STEM" & GROUP == "CTRL" & LOCATION == "sigma"~ 1,
            cell_type == "STEM" & status == "not involved" & GROUP == "CD" & LOCATION == "sigma" ~ -1,
            TRUE ~ 0),
        c28_diff_colon_H_CD = case_when( # Number 28
            cell_type == "DIFF" & status == "healthy" & GROUP == "CTRL" & LOCATION == "sigma"~ 1,
            cell_type == "DIFF" & status == "healthy" & GROUP == "CD" & LOCATION == "sigma" ~ -1,
            TRUE ~ 0)
    ) %>%
    dplyr::select(-colname, -cell_type, -status, -AaD, -age, -GROUP, -LOCATION,
                  -SAMPLE, -reanalyzed, -SAMPLE2, -TYPE) %>%
    as.matrix()

design2 <- design2[, colSums(design2 != 0) > 3]
keep <- apply(design2, 2, function(x){length(unique(x)) != 2 }) # Check that there aren't empty groups
design2 <- design2[, keep]

pick <- design %*% contr.matrix # To know which samples are taken
pick <- pick[, colSums(pick != 0) > 3 ]

pick <- cbind(pick, design2)

p <- apply(pick, 2, table)
p[lengths(p) == 2] <- lapply(p[lengths(p) == 2], function(x){
    y <- c(x, "0" = 0)
    y[c("-1", "0", "1")]})
samples_used <- t(simplify2array(p))
write.csv(samples_used[, c(3, 1)], "processed/samples_used_new.csv")

d2 <- pick
d2[d2 == 1] <- 2
d2[d2 == -1] <- 1
d2[d2 == 0] <- NA
write.csv(d2, file = "processed/comparisons.csv", row.names = TRUE)

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



# Evaluation ####

vfit <- lmFit(voom_rna, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)


dt <- decideTests(efit, lfc = log2(1.5), adjust.method = "fdr")
dtt <- summary(dt)
dtt

# Prepare for GET ####
# perl ~/Documents/projects/GETS/gets.pl --matrix=./processed/STEM_new.tsv --geneinfo=./processed/gene_info_new.tsv --sampleinfo=./processed/STEM_samples_info_new.tsv --colors=./processed/colors_info_new.tsv --output=./processed/STEM_GETS_new --center=TRUE --overwrite=TRUE
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


write.table(gene_info2[, -1], "processed/gene_info_new.tsv", sep = "\t", row.names = FALSE,
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
write.table(expr2_csv, file = "processed/new.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
stem <- cbind.data.frame(genes = rownames(expr_csv), expr_csv[, meta_o$cell_type == "STEM"])
write.table(stem,
            file = "processed/STEM_new.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")
diff <- cbind.data.frame(genes = rownames(expr_csv), expr_csv[, meta_o$cell_type == "DIFF"])
write.table(diff,
            file = "processed/DIFF_new.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")

meta_o %>%
    dplyr::filter(colname %in% colnames(expr2_csv)) %>%
    dplyr::select(-colname) %>%
    t() %>%
    write.table(file = "processed/samples_info_new.tsv", sep = "\t",
                row.names = TRUE, quote = FALSE, col.names = FALSE, na = "")
meta_o %>%
    dplyr::filter(colname %in% colnames(stem)[-1]) %>%
    dplyr::select(-colname) %>%
    t() %>%
    write.table(file = "processed/STEM_samples_info_new.tsv", sep = "\t",
                row.names = TRUE, quote = FALSE, col.names = FALSE, na = "")
meta_o %>%
    dplyr::filter(colname %in% colnames(diff)[-1]) %>%
    dplyr::select(-colname) %>%
    t() %>%
    write.table(file = "processed/DIFF_samples_info_new.tsv", sep = "\t",
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
    write.table(file = "processed/colors_info_new.tsv", quote = FALSE, sep = "\t",
                col.names = FALSE, row.names = FALSE)

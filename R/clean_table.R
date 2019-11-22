library("here")
library("edgeR")
library("tidyr")
library("ggplot2")
library("ggrepel")
library("dplyr")
library("org.Hs.eg.db")
library("raster") # For the cv
library("sva")

counts <- read.table(here("data", "ORGANOIDS_STAR_RSEM_GENES.txt"), check.names = FALSE)
isos <- read.table(here("data", "ORGANOIDS_STAR_RSEM_ISOFORMS.txt"), check.names = FALSE)

# Extract the batch
batch <- gsub(".*_", "", gsub("\\..*", "", colnames(counts)))
# Barplot to control the amount of mapped counts per sample at different batch
barplot(sort(colSums(counts)), col = as.factor(batch), border = NA)


samples <- gsub("_.*", "", colnames(counts))
uniq_samples <- unique(samples)
genes <- gsub("\\..*", "", rownames(counts))


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
umat2 <- group_samples(isos, uniq_samples)
umat2 <- umat2[rowSums(umat2) != 0, ]

# Barplot to control the amount of mapped counts per sample
barplot(sort(colSums(umat)), border = "grey")

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

meta <- readxl::read_excel("data/variables CD cohort.xlsx") %>%
    dplyr::select(rename_cols) %>%
    filter(colname %in% uniq_samples & GROUP != "UC") %>%
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

expr <- umat[, meta$colname]
expr2 <- umat2[, meta$colname]

# Checking if renalized sequence are better in number of aligned sequences.
counts_expr <- colSums(expr)
barplot(counts_expr, border = "grey", col = meta$reanalyzed)
#Scaled to the mean
barplot(scale(counts_expr)[, 1], border = "grey", col =meta$reanalyzed)

meta %>%
    subset(GROUP != "CTRL") %>%
    group_by(TYPE) %>%
    summarise(min = min(age, na.rm =TRUE),
              max = max(age, na.rm =TRUE))

subset(meta, age == 18)

meta %>%
    subset(GROUP != "CTRL") %>%
    mutate(age = as.numeric(age), AaD = as.numeric(AaD)) %>%
    ggplot()+
    geom_abline(slope = 1, intercept = 0) +
    geom_point(aes(AaD, age, shape = TYPE, col = TYPE)) +
    theme_minimal()

remove_coefficient <- meta %>%
  dplyr::group_by(GROUP, cell_type, TYPE, LOCATION) %>%
  dplyr::count(sort = TRUE) %>%
  dplyr::filter(n <= 2) %>%
  dplyr::select(cell_type, GROUP, TYPE, LOCATION) %>%
  unique() %>%
  paste(collapse = "_")


# Remove genes from Y
entrez_chr <- as.data.frame(org.Hs.egCHR)
entrez_chr <- subset(entrez_chr, chromosome %in% "Y")
ens <- mapIds(org.Hs.eg.db, keys = entrez_chr$gene_id, keytype = "ENTREZID", column = "ENSEMBL")
ens_unique <- unique(ens)
genes_entrez <- gsub("\\..*", "", rownames(expr))
# Expression without genes on Y
expr_asexual <- expr[-which(genes_entrez %in% ens_unique), ]


ens <- mapIds(org.Hs.eg.db, keys = entrez_chr$gene_id, keytype = "ENTREZID", column = "ENSEMBLTRANS")
ens_unique <- unique(ens)
genes_ensembl <- gsub("\\..*", "", rownames(expr2))
# Expression without genes on Y
expr2_asexual <- expr2[-which(genes_ensembl %in% ens_unique), ]

old_sa <- expr2_asexual[, endsWith(colnames(expr2_asexual), "V")]
new_sa <- expr2_asexual[, !endsWith(colnames(expr2_asexual), "V")]

y_old <- DGEList(old_sa)
y_new <- DGEList(new_sa)
y <- DGEList(expr2_asexual)
# filtering
# Filtrar també que no siguin els del gen X e Y



keep <- filterByExpr(y)
keep_old <- filterByExpr(y_old)
keep_new <- filterByExpr(y_new)
# stopifnot(sum(keep) == 13958) # With all the samples Salmon
# stopifnot(sum(keep) == 13973) # With subset of only relevant samples Salmon
# stopifnot(sum(keep) == 14030) # With subset of only relevant samples STAR
# stopifnot(sum(keep) == 13996) # With subset of only relevant samples STAR without Y genes
stopifnot(sum(keep) == 35608) # By isoform without Y genes
stopifnot(sum(keep_new) == 36277)
stopifnot(sum(keep_old) == 33830)

y <- y[keep, ]
y_new <- y_new[keep_new, ]
y_old <- y_old[keep_new, ]
y_norm <- calcNormFactors(y)
y_norm_new <- calcNormFactors(y_new)
y_norm_old <- calcNormFactors(y_old)


y_voom <- voom(y_norm, plot = TRUE, normalize.method = "quantile") # v$E <-normalised matrix
y_voom_new <- voom(y_norm_new, plot = TRUE, normalize.method = "quantile") # v$E <-normalised matrix
y_voom_old <- voom(y_norm_old, plot = TRUE, normalize.method = "quantile") # v$E <-normalised matrix

# PCAs
pdf("plots/PCAs_filtered_isoforms.pdf")
plotPCA(t(y_voom$E), meta$reanalyzed)
plotPCA(t(y_voom$E)[meta$cell_type == "STEM", ], meta$reanalyzed[meta$cell_type == "STEM"])
plotPCA(t(y_voom$E)[meta$cell_type == "STEM", ], meta$LOCATION[meta$cell_type == "STEM"])
plotPCA(t(y_voom$E), meta$cell_type)
plotPCA(t(y_voom$E), meta$TYPE)
plotPCA(t(y_voom$E), meta$SAMPLE)
plotPCA(t(y_voom$E), meta$LOCATION)
plotPCA(t(y_voom$E), paste(meta$reanalyzed, meta$LOCATION, sep = "_"))
plotPCA(t(y_voom$E), meta$GROUP)
plotPCA(t(y_voom$E), paste(meta$reanalyzed, meta$TYPE, sep = "_"))
plotPCA(t(y_voom$E), paste(meta$cell_type, meta$TYPE, sep = "_"))
dev.off()

# Design
grouping <- apply(meta[, c("cell_type", "GROUP", "TYPE", "LOCATION")], 1,
                paste0, collapse = "_")
grouping <- as.factor(grouping)
design <- model.matrix(~0 + grouping, grouping)
colnames(design) <- levels(grouping)
# design <- cbind(design, Reanalyzed = as.numeric(meta$reanalyzed))
samples_d <- model.matrix(~0 + SAMPLE2, data = meta[, "SAMPLE2", drop = FALSE])
design <- cbind(design, samples_d)

# Helper functions for contrasts
paste_plus <- function(...){paste(..., collapse = " + ")}
paste_ab <- function(a, b){paste(a, "- (", b, ")")}
grep_level <- function(x){
    # Force to actually look up on the design being used!
    # This makes it dependent on the order it is run!!
    a <- force(colnames(design))
    grep(x, a, value = TRUE)
}
contrasting <- function(x, y) {
    paste_ab(paste_plus(grep_level(x)),
             paste_plus(grep_level(y)))
}

contrasts <- c(
    diff_stem = contrasting("^STEM", "^DIFF"), # Siempre por separado
    CD_CTRL = contrasting("_CD_", "_CTRL_"),
    STEM__ileum_sigma = contrasting("STEM.*ileum", "STEM.*sigma"), # STEM ileum vs STEM sigma
    DIFF__ileum_sigma = contrasting("DIFF.*ileum", "DIFF.*sigma"),# DIFF ileum vs DIFF sigma
    pediatric_adult = contrasting("pediatric", "adult"),
    STEM__pediatric_adult = contrasting("STEM.*pediatric", "STEM.*adult"),
    DIFF__pediatric_adult = contrasting("DIFF.*pediatric", "DIFF.*adult"),
    ileum_sigma = contrasting("ileum", "sigma"),
    ileum__STEM_DIFF = contrasting("STEM.*ileum", "DIFF.*ileum"), # OK
    sigma__STEM_DIFF = contrasting("STEM.*sigma", "DIFF.*sigma"), # OK
    sigma_DIFF__pediatric_adult = contrasting("DIFF.*pediatric_sigma",
                                              "DIFF.*adult_sigma"),
    sigma_STEM__pediatric_adult = contrasting("STEM.*pediatric_sigma",
                                              "STEM.*adult_sigma"),
    sigma__adult_pediatric = contrasting("adult_sigma", "pediatric_sigma"),
    ileum__adult_pediatric = contrasting("adult_ileum", "pediatric_ileum"),
    sigma_diff__pediatric_adult = contrasting("DIFF.*pediatric_sigma",
                                              "DIFF.*adult_sigma"),

    DIFF_ileum_CD__pediatric_adult = contrasting("DIFF_CD_pediatric_ileum",
                                                 "DIFF_CD_adult_ileum"),
    STEM_ileum_CD__pediatric_adult = contrasting("STEM_CD_pediatric_ileum",
                                                 "STEM_CD_adult_ileum"),
    STEM_ileum__pediatric_adult = contrasting("STEM_.*_pediatric_ileum",
                                                 "STEM_.*_adult_ileum"),
    DIFF_ileum__pediatric_adult = contrasting("DIFF_.*_pediatric_ileum",
                                                 "DIFF_.*adult_ileum"),
    STEM_sigma__pediatric_adult = contrasting("STEM_.*_pediatric_sigma",
                                                 "STEM_.*_adult_sigma"),
    DIFF_sigma__pediatric_adult = contrasting("DIFF_.*_pediatric_sigma",
                                                 "DIFF_.*_adult_sigma")
)

contr.matrix <- makeContrasts(contrasts = contrasts, levels = colnames(design))
# contr.matrix
colnames(contr.matrix) <- names(contrasts)


comb2 <- ComBat(y_voom$E, meta$reanalyzed, design[, -c(ncol(design), 2)])
comb2 <- ComBat(y_voom$E, meta$reanalyzed)

pdf("plots/ComBat2.pdf")
plotPCA(t(comb2), meta$reanalyzed)
# plotPCA(t(comb2)[meta$cell_type == "STEM", ], meta$reanalyzed[meta$cell_type == "STEM"])
# plotPCA(t(comb2)[meta$cell_type == "STEM", ], meta$LOCATION[meta$cell_type == "STEM"])
plotPCA(t(comb2), meta$cell_type)
plotPCA(t(comb2), meta$TYPE)
plotPCA(t(comb2), meta$SAMPLE)
plotPCA(t(comb2), meta$LOCATION)
plotPCA(t(comb2), paste(meta$reanalyzed, meta$LOCATION, sep = "_"))
plotPCA(t(comb2), meta$GROUP)
plotPCA(t(comb2), paste(meta$reanalyzed, meta$TYPE, sep = "_"))
plotPCA(t(comb2), paste(meta$cell_type, meta$TYPE, sep = "_"))
dev.off()

vfit <- lmFit(y_voom, design, ndups = 0)

remove <- apply(design[!endsWith(meta$colname, "V"), ], 2, function(x){length(unique(x))})

vfit_new <- lmFit(y_voom_new, design[!endsWith(meta$colname, "V"), remove > 1], ndups = 0)
vfit_old <- lmFit(y_voom_old, design[endsWith(meta$colname, "V"), remove > 1], ndups = 0)
# Contrasts are here!!
vfit_new <- contrasts.fit(vfit_new, contrasts = contr.matrix[colnames(design)[remove > 1], ])
vfit_old <- contrasts.fit(vfit_old, contrasts = contr.matrix[colnames(design)[remove > 1], ])
efit <- eBayes(vfit)
efit_new <- eBayes(vfit_new)
efit_old <- eBayes(vfit_old)
plotSA(efit, main="Final model: Mean-variance trend") # We can see a plot with slope 0
plotSA(efit_new, main="Final model: Mean-variance trend") # We can see a plot with slope 0
plotSA(efit_old, main="Final model: Mean-variance trend") # We can see a plot with slope 0

#  To verify the contrasts is:
#  efit$contrasts
stopifnot(sum(efit$contrasts[, "diff_stem"]) == 0)
rownames(design) <- meta$colname
res <- design %*% efit$contrasts
stopifnot(unname(table(res[, "diff_stem"])) == unname(table(meta$cell_type)))


dt <- decideTests(efit_new, lfc = log2(1.25), adjust.method = "none")
dtt <- summary(dt)
dtt

## Split the data in two for GETS ####

# Map names to HUGO
genes <- gsub("\\..*", "", rownames(dt))
genes_n <- gsub(".*\\.[0-9]*_", "", rownames(dt))
gene_names <- mapIds(org.Hs.eg.db, keys = genes, keytype = "ENSEMBLTRANS", "SYMBOL")

gene_info <- data.frame(rownames = rownames(dt), genes, gene_names, genes_n)

# Adding data from the comparisons
tt_all <- lapply(seq_along(contrasts), function(x){
  tt <- topTable(efit, coef = x, number = Inf)
  colnames(tt) <- paste0(names(contrasts)[x], "_", colnames(tt))
  tt
})

tt_all_m <- do.call(cbind, tt_all)
tt_all_m <- tt_all_m[, !grepl("_B$|_AveExpr$|_t$", colnames(tt_all_m))]

m_data <- do.call(cbind, tt_all_m)
gene_info2 <- cbind(gene_info, m_data)


r <- order(meta$TYPE, meta$LOCATION, meta$GROUP, meta$cell_type)
x <- y_voom$E[, r]
cv <- apply(x, 1, cv)
sel_genes <- rank(abs(cv)) > length(cv)/2 # 50% of genes with more variation
x <- cbind(rownames = rownames(x), x)
x <- x[sel_genes, ]

gene_info2 <- gene_info2[sel_genes, ]

pvals <- grep("P.Value$", colnames(gene_info2))
signif <- apply(gene_info2[, pvals], 2, function(x){x < 0.05})
sign_df <- apply(gene_info2[, pvals-1], 2, sign)
s <- signif*sign_df
st <- apply(s, 2, function(x){
  y <- x
  y[x == -1] <- "DW"
  y[x == 1] <- "UP"
  y[x == 0] <- ""
  y})
colnames(st) <- gsub("_logFC", "_DEG", colnames(st))

gene_info3 <- cbind(gene_info2, st)
write.table(gene_info3, "processed/gene_info.tsv", sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)

meta <- meta[r, ]
# Split the expression data
# TODO rank rows according to a cluster, very slow
# hc <- hclust(dist(x[, -1])) # Cannot rank them   all it doesn't fit the computer...
stem <- y_voom$E[, r][, meta$cell_type == "STEM"]
diff <- y_voom$E[, r][, meta$cell_type != "STEM"]


diff <- cbind(genes = rownames(diff), diff)
stem <- cbind(genes = rownames(stem), stem)
write.table(diff, file = "processed/diff.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(stem, file = "processed/stem.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(x, file = "processed/all.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

meta %>%
  filter(colname  %in% colnames(x)) %>%
  dplyr::select(-colname) %>%
  t() %>%
  write.table(file = "processed/samples_info.tsv", sep = "\t",
            row.names = TRUE, quote = FALSE, col.names = FALSE)

val <- sapply(meta[2:ncol(meta)], unique)

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
  write.table(file = "processed/colors_info.tsv", quote = FALSE, sep = "\t",
              col.names = FALSE, row.names = FALSE)

###

# Other plots ####

as.data.frame(dtt) %>%
  filter(Var1 != "NotSig") %>%
  group_by(Var2) %>%
  summarise(diff_genes = sum(Freq)) %>%
  arrange(-diff_genes) %>%
  View()


genes <- list(dw_sigma_adults = c("PMS1", "NEIL3", "SFR1"), #
              up_ilum_adults = c("ADH4", "TTR", "PLBD1")) # ADH4 and TTR fails the quality


# Compare with the microarrays ####
microarrays <- read.table("data/nov072019_Fulldata_organoides_nonFiltered.txt",
                        check.names = FALSE, row.names = 1, sep = "\t")
microarray_n <- read.csv("data/nov072019_Fulldata_annotation.txt",
                         check.names = FALSE, row.names = 1, sep = "\t")
selected_genes <- c("PMS1", "NEIL3", "SFR1", "PLBD1", "S100A8", "S100A9")

cvs <- apply(microarrays, 1, cv)
gs <- names(cvs)[(cvs>30)]
rows <- subset(microarray_n, probes  %in% gs | hugo  %in% selected_genes)

samples_names <- gsub("RNA ", "", meta$SAMPLE)
shared_names <- intersect(samples_names, colnames(microarrays))

g <- mapIds(org.Hs.eg.db, keys = as.character(rows$hugo), keytype = "SYMBOL",
            column = "ENSEMBLTRANS", multiVals = "list")
g <- unique(unlist(g, use.names = FALSE))

rnaseq <- t(y_voom$E[gsub(".*_(.*)-.*", "\\1", rownames(y_voom$E)) %in% toupper(as.character(rows$hugo)),
                     samples_names %in% shared_names])

stopifnot(match(shared_names, meta$SAMPLE[samples_names %in% shared_names]) == seq_len(42))
colnames(rnaseq) <- gsub(".*_(.*)-.*", "\\1", colnames(rnaseq))
micro <- t(microarrays[rows$probes, shared_names])
colnames(micro) <- rows$hugo
crs <- cor(rnaseq, micro, method = "spearman")
diag(crs[, rownames(crs)])

crs2 <- as.data.frame(crs)
crs2$rownames <- rownames(crs)
crs_long <- pivot_longer(crs2, colnames(crs))
crs_shorter <- subset(crs_long, rownames == name & !is.na(value))

g <- c("S100A8", "S100A9") # New function to just the correlations between the isoforms of S100A8
pdf("plots/residuals_isoforms3.pdf")
s <- sapply(g, function(name) {
  p <- which(colnames(rnaseq) == name)
  sapply(p, function(i) {
    browser()
    cs <- data.frame(V1 = rnaseq[, i], V2 = micro[, g[!g %in% name]])
    ct <- cor.test(cs$V1, cs$V2)
    lms <- lm(V2~V1, data = cs)

    if(!is.finite(lms$coefficients[1]) || !is.finite(lms$coefficients[2])) {
      print(lms$coefficients)
      return(NULL)
    }

    plot(cs$V1, cs$V2, xlab = paste("RNAseq", colnames(rnaseq)[i]),
         ylab = paste("Microarray", g[!g %in% name]),
         type = "n",
         main = paste("r", round(ct$estimate, 5), ":p =",
                      formatC(ct$p.value, format = "e", digits = 4)
                      )
    )
    text(cs$V1, cs$V2, labels = rownames(rnaseq))
    abline(a = lms$coefficients[1], b = lms$coefficients[2], col ="blue")

  })
})
dev.off()

s2 <- do.call(rbind, s)

plot(apply(s2, 1, sd), rowSums(abs(s2))/112, type = "n")
text(apply(s2, 1, sd), rowSums(abs(s2))/112, labels = rownames(s2))




## Volcano plots
coef <- 23
topTable(efit, coef = coef, number = Inf) %>%
    as_tibble() %>%
    mutate(Selected = case_when(logFC > 1.5 & adj.P.Val < 0.05 ~ "UP",
                                logFC < -1.5 & adj.P.Val < 0.05 ~ "DOWN",
                                TRUE ~ "no change"),
           names = gene_names) %>%
    ggplot(aes(x = logFC, y = -log10(adj.P.Val), label = names)) +
    geom_point(aes(col = Selected)) +
    scale_color_manual(name = "Differential \n regulation",
                       values = c("DOWN" = "#058983", "no change" = "grey", "UP" = "#DEB132")) +
    # Filtering to just the meaningful labels
    geom_text_repel(data = function(x) {
        filter(x, names %in% unlist(genes, use.names = FALSE))},
        hjust = 0.5, segment.size = 0.2, nudge_y = 2) +
    labs(title = colnames(efit)[coef]) +
    theme_minimal()

# Comprovar expressió els gens:
apply(expr[grep("ADH4|TTR|S100A8", rownames(expr)), endsWith(colnames(expr), "V")], 1, median)

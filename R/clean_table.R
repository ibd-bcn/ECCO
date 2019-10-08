library("here")
library("edgeR")
library("dplyr")
library("tidyr")
library("ggplot2")

counts <- read.table(here("data", "ORGANOIDS_STAR_RSEM_GENES.txt"), check.names = FALSE)

samples <- gsub("_.*", "", colnames(counts))
uniq_samples <- unique(samples)
genes <- gsub("\\..*", "", rownames(counts))


# To summarize the samples together
umat <- matrix(nrow = nrow(counts), ncol = length(uniq_samples))
for (i in 1:length(uniq_samples)) {
    umat[, i] <- apply(counts[, which(samples == uniq_samples[i])], 1, sum)
}
colnames(umat) <- uniq_samples
rownames(umat) <- rownames(counts)

meta <- readxl::read_excel(here("data", "variables CD cohort.xlsx")) %>%
    rename(colname = `AGILENT/ARRAY CODE`,
           cell_type = `STEM/DIFF`,
           status = `SAMPLE STATUS`,
           AaD = `Age at diagnosis`) %>%
    subset(colname %in% uniq_samples & GROUP != "UC") %>%
    mutate(GROUP = toupper(GROUP),
           reanalyzed = endsWith(colname, "V"))

expr <- umat[, meta$colname]

meta %>%
    subset(GROUP != "CTRL") %>%
    group_by(TYPE) %>%
    summarise(min = min(age, na.rm =TRUE),
              max = max(age, na.rm =TRUE))
# Isa!! some pediatric are older than the adults!

meta %>%
    subset(GROUP != "CTRL") %>%
    mutate(age = as.numeric(age), AaD = as.numeric(AaD)) %>%
    subset(age < AaD)
# Isa!! Age at diagnosis is older than the age of the sample??

meta %>%
    subset(GROUP != "CTRL") %>%
    mutate(age = as.numeric(age), AaD = as.numeric(AaD)) %>%
    ggplot()+
    geom_abline(slope = 1, intercept = 0) +
    geom_point(aes(age, AaD, shape = TYPE, col = TYPE)) +
    theme_minimal()

remove_coefficient <- meta %>%
    group_by(GROUP, cell_type, TYPE, LOCATION) %>%
    count(sort = TRUE) %>%
    subset(n <= 2) %>%
    select(cell_type, GROUP, TYPE, LOCATION) %>%
    unique() %>%
    paste(collapse = "_")


y <- DGEList(expr)
# filtering
keep <- filterByExpr(y)
# stopifnot(sum(keep) == 13958) # With all the samples Salmon
# stopifnot(sum(keep) == 13973) # With subset of only relevant samples Salmon
stopifnot(sum(keep) == 14030) # With subset of only relevant samples STAR
y <- y[keep, ]
y_norm <- calcNormFactors(y)

y_voom <- voom(y_norm, plot = TRUE) # v$E <-normalised matrix

# Design
levels <- apply(meta[, c("cell_type", "GROUP", "TYPE", "LOCATION")], 1,
                paste0, collapse = "_")
u_levels <- unique(levels)
u_levels <- u_levels[!u_levels %in% remove_coefficient]
design <- matrix(0, nrow = nrow(meta), ncol = length(u_levels))
colnames(design) <- u_levels

for (l in colnames(design)) {
    design[levels == l, l] <- 1
}

design <- cbind(Intercept = 1, design)
# Helper functions
paste_plus <- function(...){paste(..., collapse = " + ")}
paste_ab <- function(a, b){paste(a, "- (", b, ")")}
grep_level <- function(x){
    force(u_levels)
    grep(x, u_levels, value = TRUE)
}
contrasting <- function(x, y) {
    paste_ab(paste_plus(grep_level(x)),
             paste_plus(grep_level(y)))
}

contrasts <- c(
    diff_stem = contrasting("^STEM", "^DIFF"),
    CD_CTRL = contrasting("_CD_", "_CTRL_"),
    pediatric_adult = contrasting("pediatric", "adult"),
    STEM__pediatric_adult = contrasting("STEM.*pediatric", "STEM.*adult"),
    DIFF__pediatric_adult = contrasting("DIFF.*pediatric", "DIFF.*adult"),
    ileum_sigma = contrasting("ileum", "sigma"),
    ileum__STEM_DIFF = contrasting("STEM.*ileum", "DIFF.*ileum"),
    sigma__STEM_DIFF = contrasting("STEM.*sigma", "DIFF.*sigma"),
    ileum__STEM_DIFF = contrasting("STEM.*ileum", "DIFF.*ileum"),
    sigma__STEM_DIFF = contrasting("STEM.*sigma", "DIFF.*sigma"),
    sigma_DIFF__pediatric_adult = contrasting("DIFF.*pediatric_sigma",
                                              "DIFF.*adult_sigma"),
    sigma_STEM__pediatric_adult = contrasting("STEM.*pediatric_sigma",
                                              "STEM.*adult_sigma"),
    DIFF_ileum_vs_DIFF_sigma = contrasting("DIFF.*_ileum",
                                              "DIFF.*_sigma"),
    STEM_ileum_vs_STEM_sigma = contrasting("STEM.*_ileum",
                                              "STEM.*_sigma")
)

contr.matrix <- makeContrasts(contrasts = contrasts, levels = colnames(design))
# contr.matrix
colnames(contr.matrix) <- names(contrasts)


vfit <- lmFit(y_voom, design, ndups = 0)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend") # We can see a plot with slope 0


dt <- decideTests(efit, lfc = 1.5)
summary(dt)
# Map names to HUGO
library("org.Hs.eg.db")
gene_names <- mapIds(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", "SYMBOL")

sum(dt[, 1] != 0 & dt[, 2] != 0)
volcanoplot(efit, 1, highlight = 15, names = gene_names)

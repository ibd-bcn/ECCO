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

# Subset to just old samples
expr <- umat[, colnames(umat) %in% meta$colname[meta$reanalyzed]]
meta <- meta[meta$colname %in% colnames(expr), ]

dge_rna <- DGEList(expr)
dge_rna <- calcNormFactors(dge_rna)
voom_rna <- voom(dge_rna, plot = FALSE, normalize.method = "quantile")

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


vfit <- lmFit(voom_rna, design, ndups = 0)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)


dt <- decideTests(efit_new, lfc = log2(1.25), adjust.method = "none")
dtt <- summary(dt)
dtt

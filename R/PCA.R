library("edgeR")
library("ggplot2")
library("dplyr")
library("here")

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
    dplyr::rename(rename_cols) %>%
    filter(!GROUP %in% "UC") %>%
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

counts <- round(counts, 0)
samples <- gsub("_.*", "", colnames(counts))
uniq_samples <- unique(samples)

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
umat0 <- group_samples(counts, uniq_samples)
umat <- umat0[, colnames(umat0) %in% meta$colname]
umat <- umat[rowSums(umat) != 0, ]
meta <- meta[match(colnames(umat), meta$colname), ]

y <- DGEList(umat)
# filtering
keep <- filterByExpr(y)
y <- y[keep, ]
y <- calcNormFactors(y)

v <- voom(y) # v$E <-normalised matrix

data.class <- ifelse(endsWith(colnames(umat), "V"), "Old", "New")
data.pca <- prcomp(t(v$E), scale. = TRUE)
x <- apply(data.pca$x, 2, function(x){x/10})
df <- data.frame(x, Matrigel = data.class)

ggplot(df, aes(PC1, PC2, shape = Matrigel, col = Matrigel)) +
    geom_point(size = 3) +
    # stat_ellipse(clip = "off") +
    theme_minimal() +
    labs(x = "PC1 (38.9% explained var.)",
         y = "PC2 (11.1% explained var.)",
         title = "PCA of organoids samples")
ggsave("plots/PCA_matrigel.png")

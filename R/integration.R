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

shared_samples <- intersect(colnames(umat), meta$colname)

expr <- umat[, shared_samples]
meta <- meta[match(shared_samples, meta$colname), ]

# Microbiome ####
names_microb <- read.csv("data/FINAL SAMPLES on plate.csv")
colnames(names_microb) <- c("short_code", "code", "concentr")
microb <- read.csv("data/Partek_ECCO_Kraken_Classified_genus.txt",
                   sep = "\t", row.names = 1, check.names = FALSE)
microb_names <- t(simplify2array(strsplit(colnames(microb), "_")))
microb <- as.matrix(microb)

s <- t(simplify2array(strsplit(colnames(microb), "_")))
m <- mutate(names_microb, Sample = gsub(" BACT", "", code))
m2 <- mutate(meta,
             Sample = gsub(" ?(STEM|DIFF)", "", SAMPLE),
             Sample = gsub("RNA ", "", Sample))

Location <- model.matrix(~ 0 + LOCATION, meta)[, -1, drop = FALSE]
colnames(Location) <- "ileum"
Demographics <- stats::model.matrix(~ 0 + cell_type + GROUP, data = meta)
Demographics <- Demographics[, -2]
colnames(Demographics) <- c("DIFF", "CTRL")

v <- apply(t(expr), 2, var)
A <- list(RNAseq = t(expr)[, v != 0], Micro = t(microb[, 1:149]),
          Location = Location, Demographics = Demographics)

library("RGCCA")
A_scale <- lapply(A, scale2)
taus <- numeric(4)
taus[1] <- 0.154241341909016
taus[2] <- 0.486741614668218 # ~ Aprox! Recaclucate
taus[3:4] <- 1
library("inteRmodel")
search_model()

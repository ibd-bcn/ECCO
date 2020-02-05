
# will break!! see the pegs_organoides_2017 project (not on github though)
microarray <- read.table("/home/lrevilla/Documents/users/nuria_planell/np/users/isabella_dotti/pegs_organoides_2017/processed/gene_info.tsv", sep = "\t",
                na.strings = "", quote = "", stringsAsFactors = FALSE)
colnames(microarray) <- microarray[1, ]
microarray <- microarray[-1, ]


rnaseq <- read.table("processed/genes_juanjo.tsv", sep = "\t",
                na.strings = "", quote = "", stringsAsFactors = FALSE)
colnames(rnaseq) <- rnaseq[1, ]
rnaseq <- rnaseq[-1, ]

genes <- intersect(rnaseq$Symbol, microarray$hugo)
microarray <- microarray[microarray$hugo %in% genes, ]
rnaseq <- rnaseq[rnaseq$Symbol %in% genes, ]

microarray <- microarray[match(genes, microarray$hugo), ]
rnaseq <- rnaseq[match(genes, rnaseq$Symbol), ]

microarray <- microarray[, startsWith(colnames(microarray), "fc")]
rnaseq <- rnaseq[, startsWith(colnames(rnaseq), "fc")]


pdf("processed/correlations_microarray_rnaseq.pdf")
title <- gsub("fc_contrast[0-9]+_", "", colnames(rnaseq))
for (i in seq_len(ncol(rnaseq))) {
    m <- as.numeric(microarray[, i])
    r <- as.numeric(rnaseq[, i])
    m <- log2(abs(m))*sign(m)
    r <- log2(abs(r))*sign(r)
    if (all(is.na(m)) | all(is.na(r))) {
        next
    }

    # Skip those we can't compare
    if (i %in% c(4, 8, 13:16, 19:20, 22:23, 25:28, 32:36)) {
        next
    }
    cors <- cor.test(r, m)
    plot(r, m,
         xlab = "RNAseq", ylab = "Microarray",
         main = paste("Contrast", i, "\n", title[i]), pch = 16,
         sub = paste("Correlation:", round(cors$estimate, 2),
                     "p-value:", cors$p.value))
}
dev.off()

library("dfrtopics") # Github: agoldst/dfrtopics
library("ggplot2")
library("dplyr")
library("scales")
library("decontam")

microb <- read.csv("~/Downloads/Partek_ECCO_Kraken_Classified_phylum.txt",
                   sep = "\t", row.names = 1)
colnames(microb) <- gsub("X[0-9]+_", "", colnames(microb))
microb <- as.matrix(microb)

prop_microb <- prop.table(microb, 2)
non_empty <- rowSums(microb)
m <- dfrtopics::gather_matrix(microb)

df <- as.data.frame(t(prop_microb[c(5, 20, 26), ]))
ggplot(df) +
    geom_point(aes(Bacteroidetes, Firmicutes, size = Proteobacteria))

isC <- decontam::isContaminant(t(microb), conc = rep(500, ncol(microb)),  method = "frequency")


m %>%
    filter(value != 0) %>%
    group_by(col_key) %>%
    mutate(v = value/sum(value)*100) %>%
    ggplot() +
    geom_col(aes(col_key, v, fill = row_key)) +
    scale_y_continuous(expand = c(0, 0)) +
    # scale_fill_viridis_d() +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    labs(x = "Samples", y = "%", fill = "Phylum")
ggsave("plots/percentage_microbiome.png")

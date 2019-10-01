library(tximport)
library(gdata)
library(ggbiplot)
library(limma)
library(edgeR)
library(WriteXLS)

assigncol <- function(pheno) {
  colx <- c("blue", "red", "green", "orange", "yellow", "gray", "pink", "white", "brown")
  colxm <- matrix("gray", nrow = length(pheno))
  g <- unique(pheno)

  for (i in 1:length(g)) {
    ss <- which(pheno == g[i])
    colxm[ss] <- colx[i]
  }
  if (length(colx) > 9) {
    cat("N>9:D'ont Run")
  }

  print(data.frame(group = 1:9, color = colx))

  return(colxm)
}





plotPCA <- function(xmat, gg) {
  data.class <- gg
  data.pca <- prcomp(xmat, scale. = TRUE)
  colx <- assigncol(data.class)
  g <- ggbiplot(data.pca,
    obs.scale = 0.25, var.scale = 0.5, labels.size = 3,
    groups = data.class, var.axes = FALSE, ellipse = TRUE, circle = FALSE
  )
  # g <- g + geom_point(aes(size=3,groups=data.class))
  # g <- g + scale_color_discrete(name = '')
  # g <- g + theme(legend.direction = 'horizontal',
  #              legend.position = 'top')
  # g<- g + geom_path(data=data.pca, aes(x=x, y=y,colour=colx), size=1, linetype=1)



  g <- g + theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  )


  print(g)

  g <- ggbiplot(data.pca,
    obs.scale = 0.25, var.scale = 0.5, labels.size = 3,
    groups = data.class, var.axes = FALSE, labels = rownames(xmat), ellipse = TRUE, circle = FALSE
  )
  # g <- g + geom_point(aes(size=3,groups=data.class))
  # g <- g + scale_color_discrete(name = '')
  # g <- g + theme(legend.direction = 'horizontal',
  #              legend.position = 'top')
  # g<- g + geom_path(data=data.pca, aes(x=x, y=y,colour=colx), size=1, linetype=1)



  g <- g + theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  )


  print(g)
}

extractPC123 <- function(mydata) {
  # extract from ggbiplot
  pcobj <- prcomp(mydata, scale. = TRUE)
  choices <- 1:3
  obs.scale <- 1
  var.scale <- 1
  nobs.factor <- sqrt(nrow(pcobj$x) - 1)
  d <- pcobj$sdev
  u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = "*")
  v <- pcobj$rotation
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, FUN = "*"))
  df.u <- df.u * nobs.factor
  colnames(df.u) <- c("PC1", "PC2", "PC3")
  rownames(df.u) <- rownames(mydata)
  return(df.u)
}


tx2gene <- read.table("annot_gencodev29.txt", sep = "\t", header = F)
genedb <- read.table("geneannot.txt", sep = "\t", header = F)

equiv <- tx2gene[, c(1, 3)]
colnames(equiv) <- c("TXNAME", "GENEID")

xls <- read.xls("lista.xls")
files <- paste(as.character(xls[, 1]), "/quant.sf", sep = "")


txi <- tximport(files, type = "salmon", tx2gene = NULL, txOut = TRUE)

transcripts <- txi$counts # for transcripts


rownames(txi$counts) <- rownames(txi$abundance) <- rownames(txi$length) <- rownames(txi$counts) <- equiv[, 1]


txi.sum <- summarizeToGene(txi, tx2gene)$counts

usample <- as.character(unique(xls[, 2]))

umat <- matrix(nrow = dim(txi.sum)[1], ncol = length(usample))
for (i in 1:length(usample)) {
  umat[, i] <- apply(txi.sum[, which(xls[, 2] == usample[i])], 1, sum)
}
colnames(umat) <- usample


cnt <- umat
rownames(umat) <- rownames(txi.sum)

y <- DGEList(cnt)
# filtering
keep <- filterByExpr(y)
y <- y[keep, ]
y <- calcNormFactors(y)

v <- voom(y) # v$E <-normalised matrix
dbclin <- read.xls("CLIN.xls") # database metadata

# check merging in two first columns

fulldb <- data.frame(usample, dbclin[match(colnames(v$E), dbclin[, 1]), ])
########################

loc <- fulldb$LOCATION
type <- fulldb$TYPE
STEM <- fulldb[, 4]
STATUS <- fulldb$SAMPLE.STATUS

pdf("PCA_plots.pdf")
plotPCA(t(v$E), STEM)
plotPCA(t(v$E), paste(STEM, "_", loc, sep = ""))
dev.off()

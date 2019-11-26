library(WriteXLS)
library(gplots)
library(limma)
library(edgeR)
comps <- read.csv("processed/comparisons.csv", header = T, row.names = 1, check.names = F)
data <- read.csv("processed/old_counts.csv", header = T, row.names = 1, check.names = F)

l2fc <- function(logratio, base = 2) {
  retval <- base^(logratio)
  retval <- ifelse(retval < 1, -1 / retval, retval)
  retval
}

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

  # print(data.frame(group=1:9,color=colx))

  return(colxm)
}

plotPCA <- function(xmat, gg) {
  library(ggbiplot)
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

cnt <- round(data, 0)
# DIFF_ileum_CD__pediatric_adult no two class variance here
# DIFF_ileum__pediatric_adult no two class variance here

comps <- comps
#### only 10 first

multilimma <- function(xdat, classes, nmethod) {
  # nmethod=
  # cyclicloess
  # none
  # quantile
  fc <- matrix(nrow = dim(xdat)[1], ncol = dim(classes)[2])
  p <- matrix(nrow = dim(xdat)[1], ncol = dim(classes)[2])
  fdr <- matrix(nrow = dim(xdat)[1], ncol = dim(classes)[2])
  t <- matrix(nrow = dim(xdat)[1], ncol = dim(classes)[2])

  colnames(fdr) <- paste("fdr_", colnames(classes), sep = "")
  colnames(fc) <- paste("fc_", colnames(classes), sep = "")
  colnames(p) <- paste("p_", colnames(classes), sep = "")
  colnames(t) <- paste("t_", colnames(classes), sep = "")

  rownames(fc) <- rownames(fdr) <- rownames(p) <- rownames(t) <- rownames(xdat)
  ns <- dim(classes)[2]
  for (i in 1:ns) {
    nona <- which(!is.na(classes[, i]))
    myclassx <- na.omit(classes[, i])
    if (sum(myclassx == 2) == 1 | sum(myclassx == 1) == 1) {
      next
    }
    myxdat <- xdat[, nona]
    ### norm step
    y <- DGEList(myxdat)
    y <- calcNormFactors(y)
    pdf(paste("QC_", colnames(classes)[i], "_", nmethod, ".pdf", sep = ""))
    v <- voom(y, normalize.method = nmethod, plot = TRUE) # v$E <-normalised matrix

    cat(i, colnames(classes)[i], "here\n")
    print(table(classes[, i]))

    norm <- v$E

    BaseMean <- apply(norm, 1, mean)
    Log2R <- apply(norm[, which(myclassx == 2)], 1, mean) - apply(norm[, which(myclassx == 1)], 1, mean)

    boxplot(norm, las = 3)
    plot(BaseMean, Log2R, main = colnames(classes)[i], cex = 0.3)
    abline(h = 0, lwd = 2, col = "blue")
    # to speed only most 1000 genes
    sel1000 <- order(apply(norm, 1, var), decreasing = T)
    plotPCA(t(norm[sel1000, ]), factor(myclassx))

    fc[, i] <- l2fc(Log2R)

    mod <- model.matrix(~ factor(myclassx, levels = c("1", "2")))
    fit1 <- lmFit(norm, mod)
    eb1 <- eBayes(fit1)
    plotSA(eb1, main = "Final model: Mean-variance trend")
    p[, i] <- eb1$p.value[, 2]
    fdr[, i] <- p.adjust(p[, i], method = "BH")
    tt <- topTable(eb1, coef = 2, number = Inf)
    tt <- tt[names(Log2R), ]
    plot(tt$AveExpr, tt$logFC)
    plot(tt$logFC, Log2R, cex = 0.3)
    dev.off()
    t[, i] <- eb1$t[, 2]

    lup <- length(which(p[, i] < 0.05 & fc[, i] > 1.5))
    ldw <- length(which(p[, i] < 0.05 & fc[, i] < -1.5))

    cat("UP=", lup, ": DW=", ldw, "\n")
  }
  return(list(fc = fc, p = p, fdr = fdr, t = t, BaseMean = as.matrix(BaseMean),
              Amean = as.matrix(eb1$Amean)))
}

results <- multilimma(cnt, comps, "cyclicloess")
saveRDS(results, "processed/juanjo_res.RDS")


r <- abs(results$fc) > 1.5 & results$p < 0.05
r2 <- r
r2[results$fc < -1.5 & r] <- "DW"
r2[results$fc > 1.5 & r] <- "UP"
r2[results$fc > 1.5 & results$fdr < 0.05] <- "UUP"
r2[results$fc < -1.5 & results$fdr < 0.05] <- "DDW"
r2[r2 %in% c("FALSE", "TRUE")] <- ""
colnames(r2) <- gsub("fc_", "sign_", colnames(r2))
contrasts_prefix <- paste0("contrast", seq_len(ncol(results$fc)))
colnames(r2) <- paste(contrasts_prefix, colnames(r2))
colnames(results$fc) <- paste(contrasts_prefix, colnames(results$fc))
colnames(results$p) <- paste(contrasts_prefix, colnames(results$p))
colnames(results$t) <- paste(contrasts_prefix, colnames(results$t))
colnames(results$fdr) <- paste(contrasts_prefix, colnames(results$fdr))

dfall <- data.frame(Gene = rownames(cnt), r2, results$fc, results$p, results$t, results$fdr)
write.table(dfall, file = "processed/genes_juanjo.tsv", sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE, na = "")
WriteXLS("dfall", "FULLRESULTS_CYCLICLOESS.xls")

results <- multilimma(cnt, comps, "quantile")
dfall <- data.frame(Gene = rownames(cnt), r2, results$fc, results$p, results$t, results$fdr)
r <- abs(results$fc) > 1.5 & results$p < 0.05
r2 <- r
r2[results$fc < -1.5 & r] <- "DW"
r2[results$fc > 1.5 & r] <- "UP"
r2[results$fc > 1.5 & results$fdr < 0.05] <- "UUP"
r2[results$fc < -1.5 & results$fdr < 0.05] <- "DDW"
r2[r2 %in% c("FALSE", "TRUE")] <- ""
colnames(r2) <- gsub("fc_", "sign_", colnames(r2))

contrasts_prefix <- paste0("contrast", seq_len(ncol(results$fc)))
colnames(r2) <- paste(contrasts_prefix, colnames(r2))
dfall <- data.frame(Gene = rownames(cnt), r2, results$fc, results$p, results$t, results$fdr)

WriteXLS("dfall", "FULLRESULTS_QUANTILE.xls")

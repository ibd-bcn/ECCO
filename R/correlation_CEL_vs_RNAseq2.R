
suppressMessages(library(affy))
suppressMessages(library(limma))
suppressMessages(library(gtools))
suppressMessages(library(sva))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library(hgu133plus2hsentrezg.db))
suppressMessages(library(hgu133plus2hsentrezgcdf))
suppressMessages(library("genefilter"))
suppressMessages(library(stringr))
suppressMessages(library("gsubfn"))

options(show.error.locations = TRUE)



########################################################
########## Functions not callable by the user ##########
########## Load these before anything else  ############
########################################################


remove.genes.with.low.expression <- function(data, q)
	{
	data.means <- as.data.frame(rowMeans(data))  # RNA-seq data: removing genes with low expression ( < 10th quantile)
	min.exp <- quantile(data.means[,1], q, na.rm=TRUE)
	#data <- data[rownames(data.means)[data.means[,1] > min.exp],]
	data <- data[data.means[,1] > min.exp,]

	if (nrow(data) == 0)
		{stop("[     ERROR]\tNo genes in data\n")}

	return(data)
	}


remove.genes.with.low.CV <- function(data, q)
	{

	cv <- apply(data,1,sd)/apply(data,1,mean)
	cut <- quantile(cv, q)
	filter <- which(cv>cut)

	data <- data[filter,]

	if (nrow(data) == 0)
		{stop("[     ERROR]\tNo genes in data\n")}

	return(data)
	}


plot.CEL.vs.RNAseq <- function(CEL.data, RNAseq.data, file, xlab, ylab)
	{
	cel <- rowMeans(CEL.data)
	rnaseq <- rowMeans(RNAseq.data)

	min <- min(cel, rnaseq)
	max <- max(cel, rnaseq)

	gene.list <- intersect(rownames(CEL.data), rownames(RNAseq.data))

	pdf(file)

	x.list <- vector(, length(gene.list))
	y.list <- vector(, length(gene.list))

	plot(1, type="n", xlim=c(min,max), ylim=c(min,max), main="CEL vs RNA-seq\nComparison of means", xlab=xlab, ylab=ylab)

	for (i in 1:length(gene.list))
		{
		gene <- gene.list[i]

		x <- mean(RNAseq.data[gene,])
		y <- mean(CEL.data[gene,])

		points(x, y, pch=20, cex=0.6, col='blue')

		x.list[i] <- x
		y.list[i] <- y
		}
	
	segments(min,min,max,max, lty="dashed")

	c <- cor.test(x.list, y.list, method="pearson")

	legend("topright", legend=paste("r =", round(c$estimate,2)), cex=1, bty = 'n')

	dev.off()
	}



########################################################
########## File paths ##################################
########################################################

file.RNAseq.expression = "/media/daguilar/BIOINFO2/TNF/RNA-seq/expression/RSEM.genes.counts.results"
celpath = '/media/daguilar/BIOINFO2/TNF/CEL/DVD_merged'
file.info = '/media/daguilar/BIOINFO2/TNF/CEL/bd_BCN_tnf_biopsies_10062016.tsv'

######################################################
########## Read info file about the samples ##########
######################################################

info <- read.csv(file.info, header=TRUE, sep='\t',na.string=c("n.a.","n.a",""))

info <- info[!is.na(info$sample_location) & !is.na(info$cel_file) & !is.na(info$fastq_TMG_ID),]  # removing all rows with n.a. in the 'sample_location', 'cel_file' or 'fastq_TMG_ID' columns

info <- info[info$week == 0 | info$week == 14, ]  # removing all rows where week != 0 or 14

if (nrow(info)==0)
	{stop("[     ERROR]\tSample information is empty\n", sep="");}

########################################
########## Read RNA-seq files ##########
########################################

RNAseq.data.all <- as.matrix(read.csv(file.RNAseq.expression, header=TRUE, sep='\t'))

if (sum(!(which(colnames(RNAseq.data.all) %in% info$fastq_TMG_ID))) > 0) # If there are samples in RNAseq.data which are not present in the 'info' file
	{stop("[     ERROR]\n");}

if (grepl("counts", file.RNAseq.expression) == TRUE)  # se eliminan los decimales si trabajamos con counts (hay decimales pq RSEM calcula "expected counts")
	{
	RNAseq.data.all <- round(RNAseq.data.all, 0)
	}

####################################
########## Read CEL files ##########
####################################

CEL.files.all <- intersect(as.character(info[!is.na(info$cel_file) & info$fastq_TMG_ID %in% colnames(RNAseq.data.all), "cel_file"]), list.files(celpath))
CEL.data.all <- ReadAffy(filenames=CEL.files.all, celfile.path=celpath, cdfname="hgu133plus2hsentrezgcdf")

#######################################################################################################################
########## We keep only the samples (colnames) present in both CEL and RNA-seq data (for comparison purposes) #########
#######################################################################################################################

RNAseq.data.all <- RNAseq.data.all[, as.character(info[info$cel_file %in% colnames(CEL.data.all), "fastq_TMG_ID"])]

if (ncol(RNAseq.data.all) != length(CEL.files.all))
	{stop("[     ERROR]\n");}


###################################
########## Normalization ##########
###################################

CEL.data.rma = rma(CEL.data.all, verbose=FALSE, normalize=TRUE)
CEL.data.all = exprs(CEL.data.rma)

RNAseq.data.all <- voom(RNAseq.data.all, normalize.method="quantile")
RNAseq.data.all <- RNAseq.data.all$E


#####################################################################################################
########## We change the sample (column) names to 'sample_id' in both CEL and RNA-seq data ########## 
#####################################################################################################

for (i in 1:ncol(CEL.data.all))
	{
	colnames(CEL.data.all)[i] <-  as.character(info[info$cel_file==colnames(CEL.data.all)[i], "Sample_id"])
	}

for (i in 1:ncol(RNAseq.data.all))
	{
	colnames(RNAseq.data.all)[i] <-  as.character(info[info$fastq_TMG_ID==colnames(RNAseq.data.all)[i], "Sample_id"])
	}


###################################################
########## Annotation data for CEL files ##########
###################################################

id <- as.character(rownames(CEL.data.all))

symb <- unlist(mget(id,hgu133plus2hsentrezgSYMBOL, ifnotfound=NA))
name <- unlist(mget(id,hgu133plus2hsentrezgGENENAME, ifnotfound=NA))
nameACC <- unlist(mget(id,hgu133plus2hsentrezgACCNUM, ifnotfound=NA))
entrez <- unlist(mget(id,hgu133plus2hsentrezgENTREZID, ifnotfound=NA))

annotation <- data.frame(probe=id,gene=symb,entrez=entrez,definition=name, stringsAsFactors=FALSE)


#############################################################################################################################################
########## Gene names are used as row names in CEL data (instead of probe names) to allow for an easy comparison with RNA-seq data ########## 
#############################################################################################################################################

CEL.data.all <- CEL.data.all[!(rownames(CEL.data.all) %in% (as.character(annotation[which(is.na(annotation$gene)),"probe"]))), ]
rownames(CEL.data.all) <- as.character(annotation[which(annotation$probe %in% rownames(CEL.data.all)),"gene"])

########################################
########## Filtering CEL data ##########
########################################

CEL.data.all <- remove.genes.with.low.expression(CEL.data.all, 0.10)
CEL.data.all <- remove.genes.with.low.CV(CEL.data.all, 0.25)

############################################
########## Filtering RNA-seq data ##########
############################################

RNAseq.data.all <- remove.genes.with.low.expression(RNAseq.data.all, 0.10)
RNAseq.data.all <- remove.genes.with.low.CV(RNAseq.data.all, 0.25)

######################################################################################################################
########## Calculate difference from median (~normalization) and save matrices for visualization with Excel ########## 
######################################################################################################################

CEL.data.median.all <- CEL.data.all-apply(CEL.data.all,1,median)

RNAseq.data.median.all <- RNAseq.data.all-apply(RNAseq.data.all,1,median)

##########################################
########## Plot CEL vs RNAS-seq ##########
##########################################


plot.CEL.vs.RNAseq(CEL.data.median.all, RNAseq.data.median.all, "CEL_vs_RNAseq.median.pdf", "mean of median-normalized expression levels (RNA-seq)", "mean of median-normalized expression levels (CEL)")

plot.CEL.vs.RNAseq(CEL.data.all, RNAseq.data.all, "CEL_vs_RNAseq.pdf", "mean of normalized expression levels (TPM; RNA-seq)", "mean of normalized expression levels (CEL)")



# code pulled from vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start

# install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("pasilla")
BiocManager::install("apeglm")

# ---------------- IMPORTS ---------------- #

# import packages
library("pasilla")
library("DESeq2")

# ---------------- DATA ---------------- #

# get data
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
print(pasAnno)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# import our data
our_cts <- as.matrix(read.csv("./data/features/deseq_cts.tsv", sep="\t",row.names="target_id"))
our_coldata <- read.csv("./data/features/coldata.csv", row.names = "Run")
our_coldata$Age <- factor(our_coldata$Age)
our_coldata$PMI <- factor(our_coldata$PMI)
our_coldata$pH <- factor(our_coldata$pH)
our_coldata$brain_region <- factor(our_coldata$brain_region)
our_coldata$Disorder <- factor(our_coldata$Disorder)

# look at the data
head(our_cts,2)
head(our_coldata, 2)

# not in same order! 
rownames(coldata) <- sub("fb", "", rownames(coldata))

# the same samples
all(rownames(coldata) %in% colnames(cts))
all(rownames(our_coldata) %in% colnames(our_cts))
# but not the same order!
all(rownames(coldata) == colnames(cts))
all(rownames(our_coldata) == colnames(our_cts))

# sort to be in the same order
cts <- cts[, rownames(coldata)]
our_cts <- our_cts[, rownames(our_coldata)]
all(rownames(coldata) == colnames(cts))
all(rownames(our_coldata) == colnames(our_cts))

# ---------------- DESeqDataSet ---------------- #

# create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds

# create our DESeqDataSet object
our_dds <- DESeqDataSetFromMatrix(countData = our_cts,
                                  colData = our_coldata,
                                  design = ~ PMI + pH + Disorder)

# set up metadata
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

# ---------------- FILTERING ---------------- #

# pre-filter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# factors in R!
dds$condition <- factor(dds$condition, levels = c("untreated","treated"))

# ---------------- DEA ---------------- #

# Differential expression analysis
?DESeq
# carries out: estimation of size factors, estimation of dispersion: neg. binomial GLM
dds <- DESeq(dds)
res <- results(dds)
res

# Log fold change
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resLFC

resOrdered <- resLFC[order(resLFC$pvalue),]
sum(resOrdered$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)

# ---------------- PAPER ---------------- #

# the paper approach:

# LRT
dds <- DESeq(dds, test="LRT", reduced=~1)
res <- results(dds)

## variance stabilizing
vsd <- vst(dds, blind=FALSE)


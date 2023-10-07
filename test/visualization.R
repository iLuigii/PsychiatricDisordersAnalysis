setwd(".")


AnCg_SZ <- read.csv("AnCg_SZ.csv")
AnCg_SZ_pvalues <- subset(AnCg_SZ, select = pvalue, subset = baseMean >= 10, drop = TRUE)
hist(AnCg_SZ_pvalues, col = c("red"))

AnCg_BP <- read.csv("AnCg_BP.csv")
AnCg_BP_pvalues <- subset(AnCg_BP, select = pvalue, subset = baseMean >= 10, drop = TRUE)
hist(AnCg_BP_pvalues, col = c("red"))

AnCg_MD <- read.csv("AnCg_MD.csv")
AnCg_MD_pvalues <- subset(AnCg_MD, select = pvalue, subset = baseMean >= 10, drop = TRUE)
hist(AnCg_MD_pvalues, col = c("red"))

DLPFC_SZ <- read.csv("DLPFC_SZ.csv")
DLPFC_SZ_pvalues <- subset(DLPFC_SZ, select = pvalue, subset = baseMean >= 10, drop = TRUE)
hist(DLPFC_SZ_pvalues, col = c("blue"))

DLPFC_BP <- read.csv("DLPFC_BP.csv")
DLPFC_BP_pvalues <- subset(DLPFC_BP, select = pvalue, subset = baseMean >= 10, drop = TRUE)
hist(DLPFC_BP_pvalues, col = c("blue"))

DLPFC_MD <- read.csv("DLPFC_MD.csv")
DLPFC_MD_pvalues <- subset(DLPFC_MD, select = pvalue, subset = baseMean >= 10, drop = TRUE)
hist(DLPFC_MD_pvalues, col = c("blue"))


nAcc_SZ <- read.csv("nAcc_SZ.csv")
nAcc_SZ_pvalues <- subset(nAcc_SZ, select = pvalue, subset = baseMean >= 10, drop = TRUE)
hist(nAcc_SZ_pvalues, col = c("green"))

nAcc_BP <- read.csv("nAcc_BP.csv")
nAcc_BP_pvalues <- subset(nAcc_BP, select = pvalue, subset = baseMean >= 10, drop = TRUE)
hist(nAcc_BP_pvalues, col = c("green"))

nAcc_MD <- read.csv("nAcc_MD.csv")
nAcc_MD_pvalues <- subset(nAcc_MD, select = pvalue, subset = baseMean >= 10, drop = TRUE)
hist(nAcc_MD_pvalues, col = c("green"))


install.packages("Hmisc")
library("Hmisc")

install.packages("corrplot")
library(corrplot)

mydata <- as.matrix(read.csv("lfc_data.csv", row.names = "X"))
mydata.cor = cor(mydata, method = c("spearman"))
mydata.rcorr = rcorr(as.matrix(twoA))
mydata.coeff = mydata.rcorr$r

corrplot(mydata.coeff)

# install.packages('VennDiagram')
library(VennDiagram)

AnCg_SZ_pc <- subset(AnCg_SZ, select = X, subset = pvalue <= 0.05, drop = TRUE)
AnCg_BP_pc <- subset(AnCg_BP, select = X, subset = pvalue <= 0.05, drop = TRUE)
AnCg_MD_pc <- subset(AnCg_MD, select = X, subset = pvalue <= 0.05, drop = TRUE)

nAnCg_SZ <- length(AnCg_MD_pc)
nAnCg_BP <- length(AnCg_BP_pc)
nAnCg_MD <- length(AnCg_MD_pc)
SZ_BP = length(intersect(AnCg_SZ_pc,AnCg_BP_pc))
SZ_MD = length(intersect(AnCg_SZ_pc,AnCg_MD_pc))
BP_MD = length(intersect(AnCg_BP_pc,AnCg_MD_pc))
three <- length(intersect(AnCg_SZ_pc,intersect(AnCg_BP_pc,AnCg_MD_pc)))
iSZ_BP = SZ_BP - three
iSZ_MD = SZ_MD - three
iBP_MD = BP_MD - three
SZ <- nAnCg_SZ
BP <- nAnCg_BP
MD <- nAnCg_MD
#Make the plot
grid.newpage()                                        # Move to new plotting page
draw.triple.venn(area1 = SZ,                          # Change color of venn diagram
                 area2 = BP,
                 area3 = MD,
                 n12 = SZ_BP,
                 n23 = BP_MD,
                 n13 = SZ_MD,
                 n123 = three,
                 col = "black",
                 category = c("SZ", "BP", "MP"),
                 fill = c("red", "blue", "green"))

heatdata <- as.matrix(read.csv("AnCg_SZ.csv", row.names = "X"))
heatmap(heatdata, scale = "column")
install.packages("pheatmap")
library("pheatmap")
pheatmap(heatdata, cutree_rows = 4)

Sys.setenv('R_MAX_VSIZE'=32000000000)
Sys.getenv('R_MAX_VSIZE')
install.packages("pheatmap")
library("pheatmap")
pheatmap(heatdata, cutree_rows = 4)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

install.packages('modeltools')
library(ComplexHeatmap)
Heatmap(heatdata, 
        name = "mtcars", #title of legend
        column_title = "Patients", row_title = "Genres",
        row_names_gp = gpar(fontsize = 7) # Text size for row names
)

install.packages("cluster")
library("cluster")
set.seed(2)
pa = pam(df, k = 3)
completerecords <- na.omit(heatdata) 
max(completerecords)
heatmap(completerecords, scale = "none")

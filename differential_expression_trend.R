########## Differential gene expression: Limma trend ########################################

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017")
load("counts.RData")

library(limma)
library(edgeR)
library(dplyr)
library(biomaRt)
library(org.Hs.eg.db) 
library(AnnotationDbi,lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.4")
library(ggrepel)
source("/data/htp/A07/RNA_Seq/common_files/functions.R")
#######################################################################################

##### K562 differential gene expression ###############################################

# make DGE object
dge <- DGEList(counts = counts, lib.size = colSums(counts), samples=anno, remove.zeros = TRUE)
dge <- edgeR::calcNormFactors(dge, method = "RLE")

# filter very little expressed genes
# find CPM value that corresponds to a read count of 10 
min.cpm <- cpm(10, mean(dge$samples$lib.size))
# retain genes with at least 1 sample above this threshold 
keep <- rowSums(cpm(dge) > as.vector(min.cpm)) >= 1
table(keep)
dge <- dge[keep, ]

# recalculate norm factors
dge <- edgeR::calcNormFactors(dge, method = "RLE")

# model matrix and contrast
design<-stats::model.matrix(~0+AraC,data=anno)

cont.matrix <- makeContrasts(
  yes_no = AraCYes - AraCNo,
  levels = design)

# empirical bayes model fitting
logCPM <- cpm(dge, log=TRUE, prior.count=3)
fit <- lmFit(logCPM,design)
fit2<-contrasts.fit(fit,contrasts = cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.1, trend = TRUE)

a<-topTable(fit2, number = Inf, coef=1, adjust.method = "BH",confint = TRUE)
tmp<-a[a$adj.P.Val<=0.05,]
#####################################################################################################################
#####################################################################################################################




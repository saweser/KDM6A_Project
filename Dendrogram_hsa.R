#### Dendrogram #############################################################################
rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis")
load("normalized_upm_hsa.RData")
source("functions.R")


library(RColorBrewer)
library(gplots)

#### Dendrogram K562
upm_hsa_k<-as.matrix(upm_hsa_k)
rowVar <-rowVars(upm_hsa_k)
mv100 <- order(rowVar, decreasing = T )[1:20] # take out the 100 most variable genes

heatmap_genes <- as.data.frame(upm_hsa_k[mv100,]) #select expression values of the top500 DE genes
heatmap_genes<-getGeneID_hsa(heatmap_genes)
rownames(heatmap_genes)<-heatmap_genes$external_gene_name
heatmap_genes$external_gene_name<-NULL
heatmap_genes<-as.matrix(heatmap_genes)
hmcol<-brewer.pal(9,"Purples") # color for heatmap

pdf("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots/dendrogram_K562_hsa.pdf")
heatmap.2(heatmap_genes, trace="none", margins = c(9,9), 
          dendrogram = "column", labRow = rownames(heatmap_genes), keysize = 1.2, symkey= FALSE, density.info="none", 
          key.xlab= "Expression (UPM)",cexCol=1, key.par=list(), key.xtickfun = NULL, 
          key.ytickfun = NULL,col = hmcol, colsep=1:ncol(heatmap_genes), rowsep=1:nrow(heatmap_genes),sepwidth=c(0.01,0.01))
dev.off()

#### Dendrogram HEK
upm_hsa_h<-as.matrix(upm_hsa_h)
rowVar <-rowVars(upm_hsa_h)
mv100 <- order(rowVar, decreasing = T )[1:20] # take out the 100 most variable genes

heatmap_genes <- as.data.frame(upm_hsa_h[mv100,]) #select expression values of the top500 DE genes
heatmap_genes<-getGeneID_hsa(heatmap_genes)
rownames(heatmap_genes)<-heatmap_genes$external_gene_name
heatmap_genes$external_gene_name<-NULL
heatmap_genes<-as.matrix(heatmap_genes)
hmcol<-brewer.pal(9,"Purples") # color for heatmap

pdf("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots/dendrogram_Hek_hsa.pdf")
heatmap.2(heatmap_genes, trace="none", margins = c(4,6), 
          dendrogram = "column", labRow = rownames(heatmap_genes), keysize = 1.2, symkey= FALSE, density.info="none", 
          key.xlab= "Expression (UPM)",cexCol=1, key.par=list(), key.xtickfun = NULL, 
          key.ytickfun = NULL,col = hmcol, colsep=1:ncol(heatmap_genes), rowsep=1:nrow(heatmap_genes),sepwidth=c(0.01,0.01))
dev.off()

#### Dendrogram MM1 and MM6
upm_hsa_m<-as.matrix(upm_hsa_m)
rowVar <-rowVars(upm_hsa_m)
mv100 <- order(rowVar, decreasing = T )[1:20] # take out the 100 most variable genes

heatmap_genes <- as.data.frame(upm_hsa_m[mv100,]) #select expression values of the top500 DE genes
heatmap_genes<-getGeneID_hsa(heatmap_genes)
rownames(heatmap_genes)<-heatmap_genes$external_gene_name
heatmap_genes$external_gene_name<-NULL
heatmap_genes<-as.matrix(heatmap_genes)
hmcol<-brewer.pal(9,"Purples") # color for heatmap

pdf("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots/dendrogram_MM1_MM6_hsa.pdf")
heatmap.2(heatmap_genes, trace="none", margins = c(4,6), 
          dendrogram = "column", labRow = rownames(heatmap_genes), keysize = 1.2, symkey= FALSE, density.info="none", 
          key.xlab= "Expression (UPM)",cexCol=1, key.par=list(), key.xtickfun = NULL, 
          key.ytickfun = NULL,col = hmcol, colsep=1:ncol(heatmap_genes), rowsep=1:nrow(heatmap_genes),sepwidth=c(0.01,0.01))
dev.off()


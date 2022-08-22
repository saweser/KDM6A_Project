########## Normalization and PCA #######################################################################
# normalize count data with edgeR
# then calculcate UMIs per million (UPM)
# PCA of normalized upm

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis")
load("counts.RData")
targets_mus<-read.table("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/anno/anno_mus.txt",sep = "\t")

source("functions.R")
library(edgeR)
library(ggplot2)
library(tidyr)
library(cowplot)
library(gridExtra)
library(reshape2)
library(dplyr)
##################################################################################################################

#### normalization: edgeR : all mouse samples #################################################################################
# calculate UPM of normalized counts
targets_mus$date<-as.character(targets_mus$date)
sums <- colSums(counts_mus)

# RLE
nf <- edgeR::calcNormFactors(counts_mus, method = "RLE")
upm_mus<-t(t(counts_mus*nf)/sums* 1e+06) # transform, because R will calculate row wise

upm_mus_plot<-as.data.frame(upm_mus)
upm_mus_plot$ensembl<-rownames(upm_mus_plot)
upm_mus_long<-upm_mus_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_mus_long,aes(x=sample_name,y=log2(upm)))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="RLE") +  
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box

# TMM
nf_2 <- edgeR::calcNormFactors(counts_mus, method = "TMM")
upm_mus_tmm<-as.data.frame(t((t(counts_mus*nf_2)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_mus_tmm_plot<-upm_mus_tmm
upm_mus_tmm_plot$ensembl<-rownames(upm_mus_tmm_plot)
upm_mus_tmm_long<-upm_mus_tmm_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_mus_tmm_long,aes(x=sample_name,y=log2(upm)))
p.box.t<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="TMM") + 
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box.t

plot_grid(p.box, p.box.t,ncol=2)
ggsave("normalization_mus.eps", plot = grid.arrange(p.box,p.box.t, ncol=2), width= 150, height= 200, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
########################################################################################################################

#### PCA ########################################################################################
PCA_line_mus <- pcaFunction(mat = upm_mus, inf = targets_mus, ngenes = 500, col = "cell_line")
PCA_state_mus <- pcaFunction(mat = upm_mus, inf = targets_mus, ngenes = 500, col = "state")
PCA_date_mus <- pcaFunction(mat = upm_mus, inf = targets_mus, ngenes = 500, col = "date")

ggsave("PCA_cell_line_mus.pdf", plot = PCA_line_mus, width= 100, height= 70, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
ggsave("PCA_state_mus.pdf", plot = PCA_state_mus, width= 100, height= 70, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
ggsave("PCA_date_mus.pdf", plot = PCA_date_mus, width= 100, height= 70, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
###########################################################################################################################

##### normalize for each cell line separately ######################################################################
# variance between cell lines still very big also it is not that clearly visible in the PCA (DE genes will be twice as much with normalization together)
### Normalization: EdgeR: 2A3 ##########################################################################################

# calculate UPM of normalized counts
counts_mus_2<-counts_mus[,grepl("2A3",colnames(counts_mus))]
targets_mus_2<-targets_mus[grepl("2A3",targets_mus$sample_name),]

sums <- colSums(counts_mus_2)

# normalization: RLE
nf <- edgeR::calcNormFactors(counts_mus_2, method = "RLE")
upm_mus_2<-as.data.frame(t((t(counts_mus_2*nf)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_mus_2_plot<-upm_mus_2
upm_mus_2_plot$ensembl<-rownames(upm_mus_2_plot)
upm_mus_2_long<-upm_mus_2_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_mus_2_long,aes(x=sample_name,y=log2(upm)))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="RLE") +  
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box

# normalization TMM
nf_2 <- edgeR::calcNormFactors(counts_mus_2, method = "TMM")
upm_mus_2_tmm<-as.data.frame(t((t(counts_mus_2*nf_2)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_mus_2_tmm_plot<-upm_mus_2_tmm
upm_mus_2_tmm_plot$ensembl<-rownames(upm_mus_2_tmm_plot)
upm_mus_2_tmm_long<-upm_mus_2_tmm_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_mus_2_tmm_long,aes(x=sample_name,y=log2(upm)))
p.box.t<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="TMM") + 
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box.t

plot_grid(p.box, p.box.t,ncol=2)
ggsave("normalization_mus_2A3.eps", plot = grid.arrange(p.box,p.box.t, ncol=2), width= 150, height= 200, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
########################################################################################################################

#### PCA 2A3 ##############################################################################################################
PCA_2A3_state_mus <- pcaFunction(mat = upm_mus_2, inf = targets_mus_2, ngenes = 500, col = "state2")+
  theme_grey()+
  theme(legend.title = element_blank())

PCA_2A3_date_mus <- pcaFunction(mat = upm_mus_2, inf = targets_mus_2, ngenes = 500, col = "date")+
  theme_grey()

ggsave("PCA_2A3_state_mus.pdf", plot = PCA_2A3_state_mus, width= 150, height= 120, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
ggsave("PCA_2A3_date_mus.pdf", plot = PCA_2A3_date_mus, width= 150, height= 120, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")


# add number of UMIs, genes detected and knockdown expression to targets file
# check in PCA 
genes<-read.table("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/number_genes_detected_mus.txt",sep = "\t")
colnames(genes)<-c("number_genes","sample_name","umis","genome")

upm_mus_2_kdm6a<-upm_mus_2[rownames(upm_mus_2)=="ENSMUSG00000037369",]
upm_mus_2_kdm6a_long<-melt(upm_mus_2_kdm6a)
colnames(upm_mus_2_kdm6a_long)<-c("sample_name","kdm6a")
targets_mus_2<-left_join(targets_mus_2, genes, by="sample_name")
targets_mus_2<-left_join(targets_mus_2, upm_mus_2_kdm6a_long, by="sample_name")
rownames(targets_mus_2)<-targets_mus_2$sample_name

PCA_umis <- pcaFunction(mat = upm_mus_2, inf = targets_mus_2, ngenes = 500, col = "umis")+
  theme_grey()
PCA_umis

PCA_genes <- pcaFunction(mat = upm_mus_2, inf = targets_mus_2, ngenes = 500, col = "number_genes")+
  theme_grey()
PCA_genes

PCA_kdm6a <- pcaFunction(mat = upm_mus_2, inf = targets_mus_2, ngenes = 500, col = "kdm6a")+
  theme_grey()
PCA_kdm6a

ggsave("PCA_2A3_umis_mus.pdf", plot = PCA_umis, width= 150, height= 120, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
ggsave("PCA_2A3_genes_mus.pdf", plot = PCA_genes, width= 150, height= 120, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
ggsave("PCA_2A3_kdm6a_mus.pdf", plot = PCA_kdm6a, width= 150, height= 120, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
#######################################################################################################################

#### Normalization: EdgeR: 4B9 ##########################################################################################
# calculate UPM of normalized counts
counts_mus_4<-counts_mus[,grepl("4B9",colnames(counts_mus))]
targets_mus_4<-targets_mus[grepl("4B9",targets_mus$sample_name),]

sums <- colSums(counts_mus_4)

# normalization: RLE
nf <- edgeR::calcNormFactors(counts_mus_4, method = "RLE")
upm_mus_4<-as.data.frame(t((t(counts_mus_4*nf)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_mus_4_plot<-upm_mus_4
upm_mus_4_plot$ensembl<-rownames(upm_mus_4_plot)
upm_mus_4_long<-upm_mus_4_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_mus_4_long,aes(x=sample_name,y=log2(upm)))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="RLE") +  
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box

# normalization TMM
nf_2 <- edgeR::calcNormFactors(counts_mus_4, method = "TMM")
upm_mus_4_tmm<-as.data.frame(t((t(counts_mus_4*nf_2)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_mus_4_tmm_plot<-upm_mus_4_tmm
upm_mus_4_tmm_plot$ensembl<-rownames(upm_mus_4_tmm_plot)
upm_mus_4_tmm_long<-upm_mus_4_tmm_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_mus_4_tmm_long,aes(x=sample_name,y=log2(upm)))
p.box.t<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="TMM") + 
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box.t

plot_grid(p.box, p.box.t,ncol=2)
ggsave("normalization_mus_4B9.eps", plot = grid.arrange(p.box,p.box.t, ncol=2), width= 150, height= 200, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
########################################################################################################################

#### PCA 4B9 ##############################################################################################################
PCA_4B9_state <- pcaFunction(mat = upm_mus_4, inf = targets_mus_4, ngenes = 500, col = "state2")+
  theme_grey()+
  theme(legend.title = element_blank())

PCA_4B9_date <- pcaFunction(mat = upm_mus_4, inf = targets_mus_4, ngenes = 500, col = "date")+
  theme_grey()

ggsave("PCA_4B9_state_mus.pdf", plot = PCA_4B9_state, width= 150, height= 120, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
ggsave("PCA_4B9_date_mus.pdf", plot = PCA_4B9_date, width= 150, height= 120, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")


# add number of UMIs, genes detected and knockdown expression to targets file
# check in PCA 
upm_mus_4_kdm6a<-upm_mus_4[rownames(upm_mus_4)=="ENSMUSG00000037369",]
upm_mus_4_kdm6a_long<-melt(upm_mus_4_kdm6a)
colnames(upm_mus_4_kdm6a_long)<-c("sample_name","kdm6a")
targets_mus_4<-left_join(targets_mus_4, genes, by="sample_name")
targets_mus_4<-left_join(targets_mus_4, upm_mus_4_kdm6a_long, by="sample_name")
rownames(targets_mus_4)<-targets_mus_4$sample_name

PCA_4B9_umis <- pcaFunction(mat = upm_mus_4, inf = targets_mus_4, ngenes = 500, col = "umis")+
  theme_grey()

PCA_4B9_genes <- pcaFunction(mat = upm_mus_4, inf = targets_mus_4, ngenes = 500, col = "number_genes")+
  theme_grey()


PCA_4B9_kdm6a <- pcaFunction(mat = upm_mus_4, inf = targets_mus_4, ngenes = 500, col = "kdm6a")+
  theme_grey()


ggsave("PCA_4B9_umis_mus.pdf", plot = PCA_4B9_umis, width= 150, height= 120, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
ggsave("PCA_4B9_genes_mus.pdf", plot = PCA_4B9_genes, width= 150, height= 120, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
ggsave("PCA_4B9_kdm6a_mus.pdf", plot = PCA_4B9_kdm6a, width= 150, height= 120, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
#######################################################################################################################


##### save objects for downstream analysis ###############################################################
save(list=c("upm_mus_2","upm_mus_4","targets_mus_2","targets_mus_4"),
     file="normalized_upm_mus.RData");
##############################################################################################################

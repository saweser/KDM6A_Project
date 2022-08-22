########## Normalization and PCA #######################################################################
# normalize count data with edgeR
# then calculcate UMIs per million (UPM)
# PCA of normalized upm

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis")
load("counts.RData")
targets_hsa<-read.table("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/anno/anno_hsa.txt",sep = "\t")
targets<-read.table("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/anno/anno.txt",sep = "\t")

source("functions.R")
library(edgeR)
library(ggplot2)
library(tidyr)
library(cowplot)
library(gridExtra)
######################################################################################################################

### Normalization: EdgeR: all human samples ##########################################################################################
# human: calculate UPM of normalized counts
targets_hsa$date<-as.character(targets_hsa$date)
sums <- colSums(counts_hsa)

# normalization: RLE
nf <- edgeR::calcNormFactors(counts_hsa, method = "RLE")
upm_hsa<-as.data.frame(t((t(counts_hsa*nf)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_hsa_plot<-upm_hsa
upm_hsa_plot$ensembl<-rownames(upm_hsa_plot)
upm_hsa_long<-upm_hsa_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_hsa_long,aes(x=sample_name,y=log2(upm)))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="RLE") +  
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box

# normalization TMM
nf_2 <- edgeR::calcNormFactors(counts_hsa, method = "TMM")
upm_hsa_tmm<-as.data.frame(t((t(counts_hsa*nf_2)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_hsa_tmm_plot<-upm_hsa_tmm
upm_hsa_tmm_plot$ensembl<-rownames(upm_hsa_tmm_plot)
upm_hsa_tmm_long<-upm_hsa_tmm_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_hsa_tmm_long,aes(x=sample_name,y=log2(upm)))
p.box.t<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="TMM") + 
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box.t

p<-ggplot(upm_hsa_tmm_long,aes(x=sample_name,y=log2(upm)))
p.box.t<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="TMM") + 
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box.t


plot_grid(p.box, p.box.t,ncol=2)
ggsave("normalization.eps", plot = grid.arrange(p.box,p.box.t, ncol=2), width= 150, height= 200, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
########################################################################################################################

#### PCA human ##############################################################################################################
# PCA human
PCA_line <- pcaFunction(mat = upm_hsa, inf = targets_hsa, ngenes = 500, col = "cell_line")
ggsave("PCA_hsa.pdf", plot = PCA_line, width= 100, height= 70, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
#######################################################################################################################

##### normalize for each cell line separately ######################################################################
# variance between cell lines to big to detect differences between 
### Normalization: EdgeR: K562 ##########################################################################################

# human: calculate UPM of normalized counts
counts_hsa_k<-counts_hsa[,grepl("K562",colnames(counts_hsa))]
targets_hsa_k<-targets_hsa[grepl("K562",targets_hsa$sample_name),]

sums <- colSums(counts_hsa_k)

# normalization: RLE
nf <- edgeR::calcNormFactors(counts_hsa_k, method = "RLE")
upm_hsa_k<-as.data.frame(t((t(counts_hsa_k*nf)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_hsa_k_plot<-upm_hsa_k
upm_hsa_k_plot$ensembl<-rownames(upm_hsa_k_plot)
upm_hsa_k_long<-upm_hsa_k_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_hsa_k_long,aes(x=sample_name,y=log2(upm)))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="RLE") +  
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box

# normalization TMM
nf_2 <- edgeR::calcNormFactors(counts_hsa_k, method = "TMM")
upm_hsa_k_tmm<-as.data.frame(t((t(counts_hsa_k*nf_2)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_hsa_k_tmm_plot<-upm_hsa_k_tmm
upm_hsa_k_tmm_plot$ensembl<-rownames(upm_hsa_k_tmm_plot)
upm_hsa_k_tmm_long<-upm_hsa_k_tmm_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_hsa_k_tmm_long,aes(x=sample_name,y=log2(upm)))
p.box.t<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="TMM") + 
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box.t

plot_grid(p.box, p.box.t,ncol=2)
ggsave("normalization_hsa_K562.eps", plot = grid.arrange(p.box,p.box.t, ncol=2), width= 150, height= 200, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
########################################################################################################################

#### PCA K562 ##############################################################################################################
PCA_K562 <- pcaFunction(mat = upm_hsa_k, inf = targets_hsa_k, ngenes = 500, col = "state2")+
  theme_grey()+
  theme(legend.title = element_blank())
PCA_K562

PCA_date <- pcaFunction(mat = upm_hsa_k, inf = targets_hsa_k, ngenes = 500, col = "date")+
  theme_grey()
PCA_date

# add number of UMIs, genes detected and knockdown expression to targets file
# check in PCA 
genes<-read.table("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/number_genes_detected_hsa.txt",sep = "\t")
colnames(genes)<-c("number_genes","sample_name","genome")
umis<-as.data.frame(colSums(counts_hsa))
colnames(umis)<-"umis"
umis$sample_name<-rownames(umis)
upm_hsa_k_kdm6a<-upm_hsa_k[rownames(upm_hsa_k)=="ENSG00000147050",]
upm_hsa_k_kdm6a_long<-melt(upm_hsa_k_kdm6a)
colnames(upm_hsa_k_kdm6a_long)<-c("sample_name","kdm6a")
targets_hsa_k<-left_join(targets_hsa_k, umis, by="sample_name")
targets_hsa_k<-left_join(targets_hsa_k, genes, by="sample_name")
targets_hsa_k<-left_join(targets_hsa_k, upm_hsa_k_kdm6a_long, by="sample_name")
rownames(targets_hsa_k)<-targets_hsa_k$sample_name

PCA_umis <- pcaFunction(mat = upm_hsa_k, inf = targets_hsa_k, ngenes = 500, col = "umis")+
  theme_grey()
PCA_umis

PCA_genes <- pcaFunction(mat = upm_hsa_k, inf = targets_hsa_k, ngenes = 500, col = "number_genes")+
  theme_grey()
PCA_genes

PCA_kdm6a <- pcaFunction(mat = upm_hsa_k, inf = targets_hsa_k, ngenes = 500, col = "kdm6a")+
  theme_grey()
PCA_kdm6a

ggsave("PCA_hsa_K562.pdf", plot = PCA_K562, width= 150, height= 120, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
#######################################################################################################################
##########################################################################################################################

### Normalization: EdgeR: Hek ##########################################################################################
# human: calculate UPM of normalized counts
counts_hsa_h<-counts_hsa[,grepl("Hek",colnames(counts_hsa))]
targets_hsa_h<-targets_hsa[grepl("Hek",targets_hsa$sample_name),]

sums <- colSums(counts_hsa_h)

# normalization: RLE
nf <- edgeR::calcNormFactors(counts_hsa_h, method = "RLE")
upm_hsa_h<-as.data.frame(t((t(counts_hsa_h*nf)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_hsa_h_plot<-upm_hsa_h
upm_hsa_h_plot$ensembl<-rownames(upm_hsa_h_plot)
upm_hsa_h_long<-upm_hsa_h_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_hsa_h_long,aes(x=sample_name,y=log2(upm)))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="RLE") +  
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box

# normalization TMM
nf_2 <- edgeR::calcNormFactors(counts_hsa_h, method = "TMM")
upm_hsa_h_tmm<-as.data.frame(t((t(counts_hsa_h*nf_2)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_hsa_h_tmm_plot<-upm_hsa_h_tmm
upm_hsa_h_tmm_plot$ensembl<-rownames(upm_hsa_h_tmm_plot)
upm_hsa_h_tmm_long<-upm_hsa_h_tmm_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_hsa_h_tmm_long,aes(x=sample_name,y=log2(upm)))
p.box.t<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="TMM") + 
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box.t

plot_grid(p.box, p.box.t,ncol=2)
ggsave("normalization_hsa_Hek.eps", plot = grid.arrange(p.box,p.box.t, ncol=2), width= 150, height= 200, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
########################################################################################################################

#### PCA Hek ##############################################################################################################
PCA_Hek <- pcaFunction(mat = upm_hsa_h, inf = targets_hsa_h, ngenes = 500, col = "state2")+
  theme_grey()+
  theme(legend.title = element_blank())
PCA_Hek

PCA_Hek_date <- pcaFunction(mat = upm_hsa_h, inf = targets_hsa_h, ngenes = 500, col = "date")+
  theme_grey()
PCA_Hek_date
ggsave("PCA_hsa_Hek.pdf", plot = PCA_Hek, width= 150, height= 120, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
#######################################################################################################################
###########################################################################################################################

### Normalization: EdgeR: MM1 and MM6 ##########################################################################################
# human: calculate UPM of normalized counts
counts_hsa_m<-counts_hsa[,grepl("MM",colnames(counts_hsa))]
targets_hsa_m<-targets_hsa[grepl("MM",targets_hsa$sample_name),]

sums <- colSums(counts_hsa_m)

# normalization: RLE
nf <- edgeR::calcNormFactors(counts_hsa_m, method = "RLE")
upm_hsa_m<-as.data.frame(t((t(counts_hsa_m*nf)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_hsa_m_plot<-upm_hsa_m
upm_hsa_m_plot$ensembl<-rownames(upm_hsa_m_plot)
upm_hsa_m_long<-upm_hsa_m_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_hsa_m_long,aes(x=sample_name,y=log2(upm)))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="RLE") +  
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box

# normalization TMM
nf_2 <- edgeR::calcNormFactors(counts_hsa_m, method = "TMM")
upm_hsa_m_tmm<-as.data.frame(t((t(counts_hsa_m*nf_2)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_hsa_m_tmm_plot<-upm_hsa_m_tmm
upm_hsa_m_tmm_plot$ensembl<-rownames(upm_hsa_m_tmm_plot)
upm_hsa_m_tmm_long<-upm_hsa_m_tmm_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_hsa_m_tmm_long,aes(x=sample_name,y=log2(upm)))
p.box.t<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="TMM") + 
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box.t

plot_grid(p.box, p.box.t,ncol=2)
ggsave("normalization_hsa_MM1_MM6.eps", plot = grid.arrange(p.box,p.box.t, ncol=2), width= 150, height= 200, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
########################################################################################################################

#### PCA MM1 MM6 ##############################################################################################################
PCA_MM1_MM6 <- pcaFunction(mat = upm_hsa_m, inf = targets_hsa_m, ngenes = 500, col = "cell_line")+
  theme_grey()+
  theme(legend.title = element_blank())
PCA_MM1_MM6

ggsave("PCA_hsa_MM1_MM6.pdf", plot = PCA_MM1_MM6, width= 150, height= 120, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
###############################################################################################################################################

##### save objects for downstream analysis ###############################################################
# save only the separately normalized ones
save(list=c("upm_hsa_k","upm_hsa_h","upm_hsa_m","targets_hsa_k","targets_hsa_h","targets_hsa_m"),
     file="normalized_upm_hsa.RData");
##############################################################################################################





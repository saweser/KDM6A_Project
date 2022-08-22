########## Normalization and PCA #######################################################################
# normalize count data with edgeR
# then calculcate UMIs per million (UPM)
# PCA of normalized upm

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017")
load("counts.RData")

source("/data/htp/A07/RNA_Seq/common_files/functions.R")
library(edgeR)
library(ggplot2)
library(tidyr)
library(cowplot)
library(gridExtra)
######################################################################################################################

### Normalization: EdgeR ##########################################################################################
# calculate UPM of normalized counts
anno$date<-as.character(anno$date)
sums <- colSums(counts)

# normalization: RLE
nf <- edgeR::calcNormFactors(counts, method = "RLE")
upm<-as.data.frame(t((t(counts*nf)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_plot<-upm
upm_plot$ensembl<-rownames(upm_plot)
upm_long<-upm_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_long,aes(x=sample_name,y=log2(upm)))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="RLE") +  
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box

# normalization TMM
nf_2 <- edgeR::calcNormFactors(counts, method = "TMM")
upm_tmm<-as.data.frame(t((t(counts*nf_2)/sums)* 1e+06)) # transform, because R will calculate row wise

upm_tmm_plot<-upm_tmm
upm_tmm_plot$ensembl<-rownames(upm_tmm_plot)
upm_tmm_long<-upm_tmm_plot %>% gather(key = sample_name, value = upm, -ensembl)

p<-ggplot(upm_tmm_long,aes(x=sample_name,y=log2(upm)))
p.box.t<-p+geom_boxplot(stat="boxplot",position="dodge")+
  labs(y="log2(upm)",title="TMM") + 
  theme_grey()+
  theme(axis.title.y = element_blank())+
  coord_flip()
p.box.t

plot_grid(p.box, p.box.t,ncol=2)
ggsave("normalization.eps", plot = grid.arrange(p.box,p.box.t, ncol=2), width= 150, height= 200, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")
########################################################################################################################

#### PCA ##############################################################################################################
PCA_AraC <- pcaFunction(mat = upm, inf = anno, ngenes = 500, col = "AraC")
PCA_AraC
PCA_treatment <- pcaFunction(mat = upm, inf = anno, ngenes = 500, col = "treatment",shape = "shRNA")
PCA_treatment

PCA_treatment <- pcaFunction(mat = upm[grepl("K562_shKDM6A_3.1",colnames(upm)),], inf = anno[grepl("K562_shKDM6A_3.1",anno$sample_name),], ngenes = 500, col = "treatment")




PCA_date <- pcaFunction(mat = upm, inf = anno, ngenes = 500, col = "date")
PCA_shRNA <- pcaFunction(mat = upm, inf = anno, ngenes = 500, col = "shRNA")+
  theme(legend.title = element_blank())
PCA_shRNA2 <- pcaFunction(mat = upm, inf = anno, ngenes = 500, col = "shRNA_2")+
  theme(legend.title = element_blank())

ggsave("PCA_AraC.pdf", plot = PCA_AraC, width= 100, height= 70, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")
# ggsave("PCA_AraC.pdf", plot = PCA_AraC, width= 100, height= 70, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")
#######################################################################################################################

##### save objects for downstream analysis ###########################################################################
save(list=c("upm","anno"), file="normalized_upm.RData");
######################################################################################################################





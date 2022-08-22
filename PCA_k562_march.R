########## Normalization and PCA #######################################################################
# normalize count data with edgeR
# then calculcate UMIs per million (UPM)
# PCA of normalized upm

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018")
load("counts_k562_march.RData")

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
ggsave("normalization.eps", plot = grid.arrange(p.box,p.box.t, ncol=2), width= 150, height= 200, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/plots")
########################################################################################################################

#### PCA ##############################################################################################################
PCA_treatment <- pcaFunction(mat = upm, inf = anno, ngenes = 500, col = "treatment")
PCA_treatment
PCA_sample <- pcaFunction(mat = upm, inf = anno, ngenes = 500, col = "sample")
PCA_sample

ggsave("PCA_treatment_k562_march.pdf", plot = PCA_treatment, width= 100, height= 70, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/plots")
ggsave("PCA_sample_k562_march.pdf", plot = PCA_sample, width= 100, height= 70, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/plots")
#######################################################################################################################





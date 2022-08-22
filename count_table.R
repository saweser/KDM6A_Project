##### count table ##############################################################################################################

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018")
#library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)

AllCounts<-readRDS("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/run_zUMIs/zUMIs_output/expression/k562_March.dgecounts.rds")
barcode<-read.table("/data/ngs/Barcode-Annotations/SCRBseq_96_setA.txt",sep = "\t",header = F)
anno<-read.table("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/anno.txt", header = T)

###############################################################################################################################

#### annotation ################################################################################################################
# join barcode and anno
colnames(barcode)<-c("T961s","well","barcode")
anno<-left_join(anno, barcode, by= "well")
anno$T961s<-NULL
rownames(anno)<-anno$name
##################################################################################################################################

##### counts: downsampled UMI counts #############################################################################################
names(AllCounts$umicount$exon$downsampling$downsampled_)
counts_dgC<-AllCounts$umicount$exon$downsampling$downsampled_
counts<-as.matrix(counts_dgC)
counts<-as.data.frame(counts)
# one sample got filtered out


#c<-AllCounts$readcount$exon$downsampling$downsampled_
#c<-as.matrix(c)
#c<-as.data.frame(c)
#a<-colSums(c)

# add sample name to counts
counts$gene<-rownames(counts)
counts.long <- reshape2::melt(counts)
colnames(counts.long)<-c("gene","barcode","value")

counts.long<-left_join(counts.long, anno,by="barcode")
counts.long$barcode<-NULL
counts.long$well<-NULL

counts<-reshape2::dcast(counts.long,gene~name,value.var = "value")
rownames(counts)<-counts$gene
counts$gene<-NULL
##########################################################################################################

##### same oder, same number of samples ##################################################################
# same number of samples
setdiff(rownames(anno),colnames(counts))
anno<-anno[!anno$name=="K562_F7_150218",]
# make sure that the rownames of annotation and the colanmes of counts are in the same oder
# check with 
colnames(counts)==rownames(anno)
# order anno
anno<-anno[order(anno$name),]

##### save objects for downstream analysis ###############################################################
save(list=c("counts","anno"),
     file="counts_k562_march.RData");
write.csv(counts, file = "/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/counts_k562_march2018.csv")
write.table(counts, file = "/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/counts_k562_march2018.txt",sep = "\t",row.names = TRUE, col.names = TRUE)
write.table(anno, file = "/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/anno.txt", col.names = TRUE, row.names = TRUE)
##############################################################################################################


counts<-read.table(file = "/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/counts_k562_march2018.txt")







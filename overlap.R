rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie")


F7_WT2_p<-read.table("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/F7_WT2_p.csv", sep = ",", header=T)
A4_WT2_p<-read.table("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/A4_WT2_p.csv", sep = ",", header=T)
C5_WT2_p<-read.table("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/C5_WT2_p.csv", sep = ",", header=T)

overlap<-read.table("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/overlap_A4WT2_F7WT2_C5WT2.csv", sep = ",", header=T)

a<-intersect(F7_WT2_p$external_gene_name, A4_WT2_p$external_gene_name)
b<-intersect(a, C5_WT2_p$external_gene_name)


overlap_a4f7<-read.table("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/overlap_A4WT2_F7WT2.csv", sep = ",", header=T)

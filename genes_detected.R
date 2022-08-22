##### genes detected and UMI count ###############################################################################
# zUMis was given the barcodes, so in stats the plots show the 48 barcodes only already
# therefore not necessary to run this script

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017")
load("counts.RData")

library(ggplot2)
library(gridExtra)
library(dplyr)
library(cowplot)

genecounts<-read.table("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/K562/zUMIs_output/stats/K562.genecounts.txt",header = TRUE)
umicounts<-read.table("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/K562/zUMIs_output/stats/K562.UMIcounts.txt",header = TRUE)
features<-read.table("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/K562/zUMIs_output/stats/K562.features.txt",header = TRUE)
###################################################################################################################

#### genes detected and umi counts ####################################################################################
colnames(genecounts)<-c("genes","barcode","featureType")
colnames(umicounts)<-c("umis","barcode","featureType")

number<-left_join(genecounts, umicounts, by=c("barcode","featureType"))
number<-left_join(anno, number, by="barcode")

# boxplot genes detected
p<-ggplot(number, aes(x= featureType, y=genes))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge",fill=c("turquoise4","skyblue4","plum4"))+
  labs(y="Genes detected") +  
  theme(axis.title.x = element_blank())
p.box

# boxplot umis
p<-ggplot(number, aes(x= featureType, y=umis))
p.box.umi<-p+geom_boxplot(stat="boxplot",position="dodge",fill=c("turquoise4","skyblue4","plum4"))+
  labs(y="UMIs") +  
  theme(axis.title.x = element_blank())
p.box.umi

plot_grid(p.box, p.box.umi,ncol=2)

ggsave("genes_det_k562.eps", plot = grid.arrange(p.box,p.box.umi, ncol=2), width= 200, height= 100, 
units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")

# table 
number2<-number
number2$well<-NULL
number2$barcode<-NULL

write.csv(number2, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/gene_umi_count_k562.csv",row.names = FALSE)
##############################################################################################################

#### readfeatures #############################################################################################
colnames(features)<-c("barcode","NumberOfReads","ReadsTotalPerCell","FractionReads","AssignmentType")
features<-left_join(anno, features, by="barcode")
features$AssignmentType <- factor(features$AssignmentType, levels = c("exon","intron","Ambiguity","Intergenic","Unmapped"))

# boxplot
p<-ggplot(features, aes(x= AssignmentType, y=FractionReads))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge",fill=c("turquoise4","skyblue4","plum4","palegreen4","slateblue4"))+
  labs(y="Fraction of reads") +  
  theme(axis.title.x = element_blank())
p.box

ggsave("read_features.box.eps", plot = p.box, width= 150, height= 70, 
       units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")


# mean of the AssignmentTypes
mean<-aggregate(features[, 11], list(features$AssignmentType), mean)
mean$Group.1 <- factor(mean$Group.1, levels = c("Unmapped","Intergenic","Ambiguity","intron","exon"))

# barplot
p<-ggplot(mean,aes(y=x, x=1, fill=Group.1)) 
p.bar<-p+geom_bar(stat="identity",position="stack")+ coord_flip()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.4,linetype = "solid"))+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")+
  scale_fill_manual(values=c("grey47","orange","yellow","seagreen4","royalblue4"),guide=guide_legend(reverse=T))
p.bar

ggsave("read_features.eps", plot = p.bar, width= 150, height=50, 
       units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")

 













########## Expression of knockdown gene ###########################################################
# in K562 KDM6A was knocked down with three shRNAs

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017")
load("normalized_upm.RData")

library(ggsignif)
library(reshape2)
###################################################################################################

###### K562: KDM6A expression ##############################################################################################
# enseml gene id for KDM6A: ENSG00000147050
upm_k<-upm[rownames(upm)=="ENSG00000147050",]
upm_k<-upm_k[,!grepl("AraC",colnames(upm_k))]
upm_k_long<-melt(upm_k)
colnames(upm_k_long)<-c("sample_name","value")

new<-left_join(upm_k_long,anno, by="sample_name")
new$shRNA_2<-relevel(as.factor(new$shRNA_2),"shRenilla")

# boxplot
p<-ggplot(new,aes(x=shRNA_2,y=value))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge",fill="turquoise4")+
  geom_point()+
  geom_signif(comparisons = list(c("shGFP","shKDM6A_3.1"),c("shGFP","shKDM6A_4"),c("shGFP","shKDM6A_7")), y_position = c(7, 7.5, 8),test = "t.test")+
  labs(y="KDM6A expression (UPM)") +  
  theme(axis.title.x = element_blank())
p.box

#### histogram
p<-ggplot(new[c(1:6),],aes(x=value))
p.hist<-p+geom_histogram(bins = 5)+
  theme(axis.title.x = element_blank())
p.hist

qqnorm(new[1:6,]$value)
qqline(new[1:6,]$value,col = 2)

qqnorm(new[7:12,]$value)
qqline(new[7:12,]$value,col = 2)

qqnorm(new[13:19,]$value)
qqline(new[13:19,]$value,col = 2)


ggsave("K562_KDM6A_expression.eps", plot = p.box, width= 100, height= 100, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")
##############################################################################################################################


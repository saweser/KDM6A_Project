########## Expression of knockdown gene ###########################################################
# in 2A3 and 4B9 KDM6A was knocked down 

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis")
load("normalized_upm_mus.RData")

library(ggsignif)
library(reshape2)
###################################################################################################

###### 2A3: KDM6A expression ##############################################################################################
# enseml gene id for KDM6A: ENSMUSG00000037369 already in targets file

# boxplot
p<-ggplot(targets_mus_2,aes(x=state2,y=kdm6a))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge",fill="turquoise4")+
  geom_point()+
  geom_signif(comparisons = list(c("KO","WT")), y_position = c(8.5),test = "wilcox.test")+
  labs(y="KDM6A expression (UPM)") +  
  theme(axis.title.x = element_blank())
p.box

p<-ggplot(targets_mus_2,aes(x=sample_name,y=kdm6a,fill=state2)) 
p.bar<-p+geom_bar(stat="identity",position="stack") +coord_flip()+
  labs(y="KDM6A expression (UPM)")+
  theme(legend.title = element_blank())+
  theme(axis.title.y = element_blank())+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c("turquoise4", "#9999CC", "steelblue"), guide=guide_legend(reverse=T))

p.bar

ggsave("2A3_KDM6A_expression_box.eps", plot = p.box, width= 100, height= 100, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
ggsave("2A3_KDM6A_expression_bar.eps", plot = p.bar, width= 100, height= 100, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
##############################################################################################################################

###### 4B9: KDM6A expression ##############################################################################################
# enseml gene id for KDM6A: ENSMUSG00000037369

# boxplot
p<-ggplot(targets_mus_4,aes(x=state2,y=kdm6a))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge",fill="turquoise4")+
  geom_point()+
  geom_signif(comparisons = list(c("KO","WT")), y_position = c(8.5),test = "wilcox.test")+
  labs(y="KDM6A expression (UPM)") +  
  theme(axis.title.x = element_blank())
p.box

p<-ggplot(targets_mus_4,aes(x=sample_name,y=kdm6a,fill=state2)) 
p.bar<-p+geom_bar(stat="identity",position="stack") +coord_flip()+
  labs(y="KDM6A expression (UPM)")+
  theme(legend.title = element_blank())+
  theme(axis.title.y = element_blank())+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c("turquoise4", "#9999CC", "steelblue"), guide=guide_legend(reverse=T))

p.bar

ggsave("4B9_KDM6A_expression_box.eps", plot = p.box, width= 100, height= 100, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
ggsave("4B9_KDM6A_expression_bar.eps", plot = p.bar, width= 100, height= 100, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
##############################################################################################################################






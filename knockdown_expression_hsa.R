########## Expression of knockdown gene ###########################################################
# in K562 and Hek293T KDM6A was knocked down with two siRNAs (siRNA_32 and siRNA34)

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis")
load("normalized_upm_hsa.RData")

library(ggsignif)
library(reshape2)
###################################################################################################

###### K562: KDM6A expression ##############################################################################################
# enseml gene id for KDM6A: ENSG00000147050
upm_hsa_k<-upm_hsa_k[rownames(upm_hsa_k)=="ENSG00000147050",]
upm_hsa_k_long<-melt(upm_hsa_k)
colnames(upm_hsa_k_long)<-c("sample_name","value")

new<-left_join(upm_hsa_k_long,targets_hsa_k, by="sample_name")

# boxplot
p<-ggplot(new,aes(x=state2,y=value))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge",fill="turquoise4")+
  geom_point()+
  geom_signif(comparisons = list(c("scr","siRNA_32"),c("scr","siRNA_34")), y_position = c(6, 5.5),test = "wilcox.test")+
  labs(y="KDM6A expression (UPM)") +  
  theme(axis.title.x = element_blank())
p.box

p<-ggplot(new,aes(x=sample_name,y=value,fill=state2)) 
p.bar<-p+geom_bar(stat="identity",position="stack") +coord_flip()+
  labs(y="KDM6A expression (UPM)")+
  theme(legend.title = element_blank())+
  theme(axis.title.y = element_blank())+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c("turquoise4", "#9999CC", "steelblue"), guide=guide_legend(reverse=T))

p.bar

ggsave("K562_KDM6A_expression_box.eps", plot = p.box, width= 100, height= 100, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
ggsave("K562_KDM6A_expression_bar.eps", plot = p.bar, width= 100, height= 100, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
##############################################################################################################################

###### Hek293T: KDM6A expression ##############################################################################################
# enseml gene id for KDM6A: ENSG00000147050
upm_hsa_h<-upm_hsa_h[rownames(upm_hsa_h)=="ENSG00000147050",]
upm_hsa_h_long<-melt(upm_hsa_h)
colnames(upm_hsa_h_long)<-c("sample_name","value")

new_h<-left_join(upm_hsa_h_long,targets_hsa_h, by="sample_name")

# boxplot
p<-ggplot(new_h,aes(x=state2,y=value))
p.box<-p+geom_boxplot(stat="boxplot",position="dodge",fill="turquoise4")+
  geom_point()+
  geom_signif(comparisons = list(c("scr","siRNA_32"),c("scr","siRNA_34")), y_position = c(6, 5.6),test = "wilcox.test")+
  labs(y="KDM6A expression (UPM)") +  
  theme(axis.title.x = element_blank())

p.box

p<-ggplot(new_h,aes(x=sample_name,y=value,fill=state2)) 
p.bar<-p+geom_bar(stat="identity",position="stack") +coord_flip()+
  labs(y="KDM6A expression (UPM)")+
  theme(legend.title = element_blank())+
  theme(axis.title.y = element_blank())+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values=c("turquoise4", "#9999CC", "steelblue"), guide=guide_legend(reverse=T))

p.bar

ggsave("Hek293T_KDM6A_expression_box.eps", plot = p.box, width= 100, height= 100, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
ggsave("Hek293T_KDM6A_expression_bar.eps", plot = p.bar, width= 100, height= 100, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
####################################################################################################




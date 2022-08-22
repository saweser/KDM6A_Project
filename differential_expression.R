########## Differential gene expression: Limma ########################################

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/KDM6A_Project/Sophie/k562_Sept2017")
load("counts.RData")

library(limma)
library(edgeR)
library(dplyr)
library(biomaRt)
library(org.Hs.eg.db) 
library(AnnotationDbi)
library(ggrepel)
library(extrafont)
loadfonts()
source("/data/htp/A07/RNA_Seq/common_files/functions.R")
#######################################################################################

##### K562 differential gene expression #############################################################################
# make DGE object
dge <- DGEList(counts = counts, lib.size = colSums(counts), samples=anno, remove.zeros = TRUE)
dge <- edgeR::calcNormFactors(dge, method = "RLE")

# filter very little expressed genes
# find CPM value that corresponds to a read count of 10 
min.cpm <- cpm(10, mean(dge$samples$lib.size))
# retain genes with at least 1 sample above this threshold 
keep <- rowSums(cpm(dge) > as.vector(min.cpm)) >= 1
table(keep)
dge <- dge[keep, ]

# recalculate norm factors
dge <- edgeR::calcNormFactors(dge, method = "RLE")
#####################################################################################################################

#### model matrix and contrast: AraC yes vs no ######################################################################################
design<-stats::model.matrix(~0+AraC,data=anno)

cont.matrix <- makeContrasts(
  yes_no = AraCYes - AraCNo,
  levels = design)

# empirical bayes model fitting
v <- voom(dge,design = design, plot = TRUE)
fit <- lmFit(v,design)
fit2<-contrasts.fit(fit,contrasts = cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.1, trend = F)

a<-topTable(fit2, number = Inf, coef=1, adjust.method = "BH",confint = TRUE)
tmp<-a[a$adj.P.Val<=0.05,]
#####################################################################################################################

#### model matrix and contrast: contrasts of sh RNAs ###############################################################
design<-stats::model.matrix(~0+treatment,data=anno)
colnames(design)<-c("GFP","KDM6A_3.1","KDM6A_3.1_AraC","KDM6A_4","KDM6A_7","KDM6A_7_AraC","Renilla","Renilla_AraC")

cont.matrix <- makeContrasts(
  sh3.1_renilla= KDM6A_3.1 - Renilla,
  sh3.1_GFP = KDM6A_3.1 - GFP,
  sh4_renilla = KDM6A_4 - Renilla,
  sh4_GFP = KDM6A_4 - GFP,
  sh7_renilla = KDM6A_7 - Renilla,
  sh7_GFP = KDM6A_7 - GFP,
  renillaAraC_renilla = Renilla_AraC - Renilla,
  k3.1AraC_k3.1 = KDM6A_3.1_AraC - KDM6A_3.1,
  k7AraC_k7 = KDM6A_7_AraC - KDM6A_7,
  k3.1AraC_renillaAraC = KDM6A_3.1_AraC - Renilla_AraC,
  k7AraC_renillaAraC =  KDM6A_7_AraC - Renilla_AraC,
  renilla_GFP = Renilla - GFP,
  levels = design)

# empirical bayes model fitting
v <- voom(dge,design = design, plot = TRUE)
fit <- lmFit(v,design)
fit2<-contrasts.fit(fit,contrasts = cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.1, trend = F)

# contrast 1: shRenilla vs shKDM6A_3.1
sh3.1_renilla<-topTable(fit2, number = Inf, coef=1, adjust.method = "BH",confint = TRUE)
sh3.1_renilla_p<-sh3.1_renilla[sh3.1_renilla$adj.P.Val<=0.05,]

# contrast 2: shGFP vs shKDM6A_3.1
sh3.1_GFP<-topTable(fit2, number = Inf, coef=2, adjust.method = "BH",confint = TRUE)
sh3.1_GFP_p<-sh3.1_GFP[sh3.1_GFP$adj.P.Val<=0.05,]

# contrast 3: shRenilla vs shKDM6A_4
sh4_renilla<-topTable(fit2, number = Inf, coef=3, adjust.method = "BH",confint = TRUE)
sh4_renilla_p<-sh4_renilla[sh4_renilla$adj.P.Val<=0.05,]

# contrast 4: shGFP vs shKDM6A_4
sh4_GFP<-topTable(fit2, number = Inf, coef=4, adjust.method = "BH",confint = TRUE)
sh4_GFP_p<-sh4_GFP[sh4_GFP$adj.P.Val<=0.05,]

# contrast 5: shRenilla vs shKDM6A_7
sh7_renilla<-topTable(fit2, number = Inf, coef=5, adjust.method = "BH",confint = TRUE)
sh7_renilla_p<-sh7_renilla[sh7_renilla$adj.P.Val<=0.05,]

# contrast 6: shGFP vs shKDM6A_7
sh7_GFP<-topTable(fit2, number = Inf, coef=6, adjust.method = "BH",confint = TRUE)
sh7_GFP_p<-sh7_GFP[sh7_GFP$adj.P.Val<=0.05,]

# contrast 7: shRenilla vs shRenilla plus AraC
renillaAraC_renilla<-topTable(fit2, number = Inf, coef=7, adjust.method = "BH",confint = TRUE)
renillaAraC_renilla_p<-renillaAraC_renilla[renillaAraC_renilla$adj.P.Val<=0.05,]

# contrast 8: shKDM6A_3.1 vs shKDM6A_3.1 plus AraC
k3.1AraC_k3.1<-topTable(fit2, number = Inf, coef=8, adjust.method = "BH",confint = TRUE)
k3.1AraC_k3.1_p<-k3.1AraC_k3.1[k3.1AraC_k3.1$adj.P.Val<=0.05,]

# contrast 9: shKDM6A_7 vs shKDM6A_7 plus AraC
k7AraC_k7<-topTable(fit2, number = Inf, coef=9, adjust.method = "BH",confint = TRUE)
k7AraC_k7_p<-k7AraC_k7[k7AraC_k7$adj.P.Val<=0.05,]

# contrast 10: shRenilla plus AraC vs shKDM6A_3.1 plus AraC
k3.1AraC_renillaAraC<-topTable(fit2, number = Inf, coef=10, adjust.method = "BH",confint = TRUE)
k3.1AraC_renillaAraC_p<-k3.1AraC_renillaAraC[k3.1AraC_renillaAraC$adj.P.Val<=0.05,]

# contrast 11: shRenilla plus AraC vs shKDM6A_7 plus AraC
k7AraC_renillaAraC<-topTable(fit2, number = Inf, coef=11, adjust.method = "BH",confint = TRUE)
k7AraC_renillaAraC_p<-k7AraC_renillaAraC[k7AraC_renillaAraC$adj.P.Val<=0.05,]

# contrast 12: shRenilla vs shGFP
renilla_GFP<-topTable(fit2, number = Inf, coef=12, adjust.method = "BH",confint = TRUE)
renilla_GFP_p<-renilla_GFP[renilla_GFP$adj.P.Val<=0.05,]
#######################################################################################################################

### to send data to Shady Awad from the hematology research unit helsinki #############################################
#write.xlsx(sh7_GFP, "/data/htp/A07/KDM6A_Project/Sophie/k562_Sept2017/tables/sh7.vs.GFP.xlsx")
#aa<-read.csv("/data/htp/A07/KDM6A_Project/Sophie/k562_Sept2017/tables/sh7_GFP_p.csv", row.names = 1)
#write.xlsx(aa, "/data/htp/A07/KDM6A_Project/Sophie/k562_Sept2017/tables/sh7_GFP_p.xlsx")
#######################################################################################################################


#### get gene names: significant ones ###################################################################################################
sh3.1_renilla_p<-getGeneID_hsa(sh3.1_renilla_p)
sh3.1_GFP_p<-getGeneID_hsa(sh3.1_GFP_p)
sh4_renilla_p<-getGeneID_hsa(sh4_renilla_p)
sh4_GFP_p<-getGeneID_hsa(sh4_GFP_p)
sh7_renilla_p<-getGeneID_hsa(sh7_renilla_p)
sh7_GFP_p<-getGeneID_hsa(sh7_GFP_p)
renillaAraC_renilla_p<-getGeneID_hsa(renillaAraC_renilla_p)
k3.1AraC_k3.1_p<-getGeneID_hsa(k3.1AraC_k3.1_p)
k7AraC_k7_p<-getGeneID_hsa(k7AraC_k7_p)
k3.1AraC_renillaAraC_p<-getGeneID_hsa(k3.1AraC_renillaAraC_p)
k7AraC_renillaAraC_p<-getGeneID_hsa(k7AraC_renillaAraC_p)
#######################################################################################################################

#### get gene names: whole table ###################################################################################################
sh3.1_renilla<-getGeneID_hsa(sh3.1_renilla)
sh3.1_GFP<-getGeneID_hsa(sh3.1_GFP)
sh4_renilla<-getGeneID_hsa(sh4_renilla)
sh4_GFP<-getGeneID_hsa(sh4_GFP)
sh7_renilla<-getGeneID_hsa(sh7_renilla)
sh7_GFP<-getGeneID_hsa(sh7_GFP)
renillaAraC_renilla<-getGeneID_hsa(renillaAraC_renilla)
k3.1AraC_k3.1<-getGeneID_hsa(k3.1AraC_k3.1)
k7AraC_k7<-getGeneID_hsa(k7AraC_k7)
k3.1AraC_renillaAraC<-getGeneID_hsa(k3.1AraC_renillaAraC)
k7AraC_renillaAraC<-getGeneID_hsa(k7AraC_renillaAraC)
#######################################################################################################################

#### decideTests ######################################################################################################
dec_tests <- decideTests(object = fit2, method = "separate", adjust.method = "BH", p.value = 0.05)
summary(dec_tests)
#######################################################################################################################

#### venn diagram #####################################################################################################
#### venn diagram of the sh RNAs with and without Arac
dec_tests_sub<-dec_tests[,grepl("AraC",colnames(dec_tests))]
dec_tests_sub<-dec_tests_sub[,1:3]

pdf("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots/venn_AraC.pdf")
limma::vennDiagram(dec_tests_sub,include = "both",circle.col = c("royalblue","seagreen2","plum","slateblue"),cex = c(1,0.8,0.8))
dev.off()

# upsetR
a<-as.data.frame(dec_tests_sub)
a[a==-1]<-1
colnames(a)<-c("shRenilla+AraC", "shKDM6A#3+AraC", "shKDM6A#7+AraC")
upset(a, sets = c("shRenilla+AraC", "shKDM6A#7+AraC"), order.by = "freq")

pdf(file = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots/upset_AraC.pdf",width=6, height=3,onefile=FALSE)
upset(a, sets = c("shRenilla+AraC", "shKDM6A#7+AraC"), order.by = "freq")
dev.off()

# overlap
overlap <- VennDiagram::calculate.overlap(list(rownames(renillaAraC_renilla[renillaAraC_renilla$adj.P.Val<=0.05,]),
                                               rownames(k3.1AraC_k3.1[k3.1AraC_k3.1$adj.P.Val<=0.05,]),
                                               rownames(k7AraC_k7[k7AraC_k7$adj.P.Val<=0.05,])))
# a1-a7: circles in venn diagram (start counting from the left to right, then second row, left to right etc.)
a5<-renillaAraC_renilla %>% dplyr::filter(rownames(renillaAraC_renilla)%in% overlap$a5)
a6<-k7AraC_k7 %>% dplyr::filter(rownames(k7AraC_k7) %in% overlap$a6)

write.csv(a5, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/overlap_AraC_shRNA_all_three.csv")
write.csv(a6, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/overlap_AraC_shRNA_3.1_7_AraC.csv")
#######################################################
#### venn diagram of shRNAs vs renilla
dec_tests_sub_2<-dec_tests[,!grepl("AraC",colnames(dec_tests))]
dec_tests_sub_2<-dec_tests_sub_2[,c(1,3,5)]

pdf("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots/venn_renilla.pdf")
limma::vennDiagram(dec_tests_sub_2,include = "both",circle.col = c("royalblue","seagreen2","plum","slateblue"),cex = c(1,0.8,0.8))
dev.off()

# upsetR
a<-as.data.frame(dec_tests_sub_2)
a[a==-1]<-1
colnames(a)<-c("shKDMA#3", "shKDM6A#4", "shKDM6A#7")
upset(a, sets = c("shKDMA#3", "shKDM6A#4", "shKDM6A#7"), order.by = "freq")


pdf(file = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots/upset_shRNA_renilla.pdf",width=6, height=3,onefile=FALSE)
upset(a, sets = c("shKDMA#3", "shKDM6A#4", "shKDM6A#7"), order.by = "freq")
dev.off()
  
  
# overlap
overlap <- VennDiagram::calculate.overlap(list(rownames(sh3.1_renilla[sh3.1_renilla$adj.P.Val<=0.05,]),
                                               rownames(sh4_renilla[sh4_renilla$adj.P.Val<=0.05,]),
                                               rownames(sh7_renilla[sh7_renilla$adj.P.Val<=0.05,])))
# a1-a7: circles in venn diagram (start counting from the left to right, then second row, left to right etc.)
a5<-sh3.1_renilla %>% dplyr::filter(rownames(sh3.1_renilla)%in% overlap$a5)
write.csv(a5, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/overlap_shRNA_renilla_all_three.csv")
######################################################
#### venn diagram of shRNAs vs GFP
dec_tests_sub_3<-dec_tests[,!grepl("AraC",colnames(dec_tests))]
dec_tests_sub_3<-dec_tests_sub_3[,c(2,4,6)]

pdf("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots/venn_GFP.pdf")
limma::vennDiagram(dec_tests_sub_3,include = "both",circle.col = c("royalblue","seagreen2","plum","slateblue"),cex = c(1,0.8,0.8))
dev.off()


# upsetR
a<-as.data.frame(dec_tests_sub_3)
a[a==-1]<-1
colnames(a)<-c("shKDM6A#3", "shKDM6A#4", "shKDM6A#7")
upset(a, sets = c("shKDM6A#3", "shKDM6A#4", "shKDM6A#7"), order.by = "freq")


pdf(file = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots/upset_shRNA_GFP.pdf",width=6, height=3,onefile=FALSE)
upset(a, sets = c("shKDM6A#3", "shKDM6A#4", "shKDM6A#7"), order.by = "freq")
dev.off()




# overlap
overlap <- VennDiagram::calculate.overlap(list(rownames(sh3.1_GFP[sh3.1_GFP$adj.P.Val<=0.05,]),
                                               rownames(sh4_GFP[sh4_GFP$adj.P.Val<=0.05,]),
                                               rownames(sh7_GFP[sh7_GFP$adj.P.Val<=0.05,])))
# a1-a7: circles in venn diagram (start counting from the left to right, then second row, left to right etc.)
a5<-sh3.1_GFP %>% dplyr::filter(rownames(sh3.1_GFP)%in% overlap$a5)
write.csv(a5, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/overlap_shRNA_GFP_all_three.csv")
#####################################################
#### vennn diagram of the shRNA 7
dec_tests_sub_4<-dec_tests[,grepl("AraC",colnames(dec_tests))]
dec_tests_sub_4<-dec_tests_sub_4[,c(1,3,5)]

pdf("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots/venn_shRNA7.pdf")
limma::vennDiagram(dec_tests_sub_4,include = "both",circle.col = c("royalblue","seagreen2","plum","slateblue"),cex = c(1,0.8,0.8))
dev.off()


pdf(file = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots/upset_shRNA_GFP.pdf",width=6, height=3,onefile=FALSE)
upset(a, sets = c("shKDMA#3", "shKDM6A#4", "shKDM6A#7"), order.by = "freq")
dev.off()



# overlap
overlap <- VennDiagram::calculate.overlap(list(rownames(renillaAraC_renilla[renillaAraC_renilla$adj.P.Val<=0.05,]),
                                               rownames(k7AraC_k7[k7AraC_k7$adj.P.Val<=0.05,]),
                                               rownames(k7AraC_renillaAraC[k7AraC_renillaAraC$adj.P.Val<=0.05,])))
# a1-a7: circles in venn diagram (start counting from the left to right, then second row, left to right etc.)
a5<-renillaAraC_renilla %>% dplyr::filter(rownames(renillaAraC_renilla)%in% overlap$a5)
write.csv(a5, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/overlap_sh7_all_three.csv")
###################################################

#### venn diagram of sh4_GFP, sh4_Renilla
dec_tests_sub_5<-dec_tests[,c(3:4)]

pdf("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots/venn_sh4_renilla_GFP.pdf")
limma::vennDiagram(dec_tests_sub_5,include = "both",circle.col = c("royalblue","seagreen2","plum","slateblue"),cex = c(1,0.8,0.8))
dev.off()

# overlap
overlap <- VennDiagram::calculate.overlap(list(rownames(sh4_renilla[sh4_renilla$adj.P.Val<=0.05,]),
                                               rownames(sh4_GFP[sh4_GFP$adj.P.Val<=0.05,])))

# with only 2 comparisons cirles are different!!!!!!!
a3<-sh4_renilla %>% dplyr::filter(rownames(sh4_renilla)%in% overlap$a3)
write.csv(a3, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/overlap_sh4_renilla_sh4_GFP.csv")
####################################################

#### venn diagram sh7_GFP und sh7_Renilla
dec_tests_sub_6<-dec_tests[,c(5:6)]

pdf("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots/venn_sh7_renilla_GFP.pdf")
limma::vennDiagram(dec_tests_sub_6,include = "both",circle.col = c("royalblue","seagreen2","plum","slateblue"),cex = c(1,0.8,0.8))
dev.off()

# overlap
overlap <- VennDiagram::calculate.overlap(list(rownames(sh7_renilla[sh7_renilla$adj.P.Val<=0.05,]),
                                               rownames(sh7_GFP[sh7_GFP$adj.P.Val<=0.05,])))

# with only 2 comparisons cirles are different!!!!!!!
a3<-sh7_renilla %>% dplyr::filter(rownames(sh7_renilla)%in% overlap$a3)
write.csv(a3, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/overlap_sh7_renilla_sh7_GFP.csv")
####################################################


##### volcano plot #######################################################################################################
# sh3.1_GFP
sh3.1_GFP<-sh3.1_GFP%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p.1 <- ggplot(data = sh3.1_GFP, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  labs(title="sh3.1_GFP")+
  geom_text_repel(data = subset(sh3.1_GFP, adj.P.Val< 0.05 |adj.P.Val < 0.05 & logFC>=1 | adj.P.Val<0.05 & logFC<=-1), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6))
  
volcano.p.1 

ggsave("volcano_sh3.1_GFP.eps", plot = volcano.p.1, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")

# sh4_GFP
sh4_GFP<-sh4_GFP%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p.2 <- ggplot(data = sh4_GFP, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  labs(title="sh4_GFP")+
  geom_text_repel(data = subset(sh4_GFP, adj.P.Val< 0.001 |adj.P.Val < 0.05 & logFC>=1 | adj.P.Val<0.05 & logFC<=-1), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)) 
volcano.p.2 

ggsave("volcano_sh4_GFP.eps", plot = volcano.p.2, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")

# sh7_GFP
sh7_GFP<-sh7_GFP%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p.3 <- ggplot(data = sh7_GFP, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  labs(title="sh7_GFP")+
  geom_text_repel(data = subset(sh7_GFP, adj.P.Val< 0.00000001 | adj.P.Val < 0.05 & logFC>=2.5 | adj.P.Val<0.05 & logFC<=-2.5 | external_gene_name =="SLC29A1"), 
                  aes(label = external_gene_name), 
                  size = 4, family="Arial",
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 16, colour = "black",family = "Arial"),
        axis.text = element_text(size = 12, colour = "black",family = "Arial"), strip.text = element_text(size = 10,family = "Arial")) 
volcano.p.3 

ggsave("volcano_sh7_GFP.pdf", plot = volcano.p.3, width= 150, height= 150, units= "mm", device = "pdf",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")

# k3.1AraC_k3.1
k3.1AraC_k3.1<-k3.1AraC_k3.1%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p.4 <- ggplot(data = k3.1AraC_k3.1, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  labs(title="k3.1AraC_k3.1")+
  geom_text_repel(data = subset(k3.1AraC_k3.1, adj.P.Val<  0.0000000001 |adj.P.Val < 0.05 & logFC>=3 | adj.P.Val<0.05 & logFC<=-3), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)) 
volcano.p.4 

ggsave("volcano_k3.1AraC_k3.1.eps", plot = volcano.p.4, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")

## k7AraC_k7
k7AraC_k7<-k7AraC_k7%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p.5 <- ggplot(data = k7AraC_k7, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  labs(title="k7AraC_k7")+
  geom_text_repel(data = subset(k7AraC_k7, adj.P.Val< 0.0000000001 |adj.P.Val < 0.05 & logFC>=3 | adj.P.Val<0.05 & logFC<=-3), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)) 
volcano.p.5 

ggsave("volcano_k7AraC_k7.eps", plot = volcano.p.5, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")


# renilla_3.1
sh3.1_renilla<-sh3.1_renilla%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p.6 <- ggplot(data = sh3.1_renilla, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  labs(title="sh3.1_renilla")+
  geom_text_repel(data = subset(sh3.1_renilla, adj.P.Val< 0.001 |adj.P.Val < 0.05 & logFC>=1 | adj.P.Val<0.05 & logFC<=-1), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)) 
volcano.p.6 

ggsave("sh3.1_renilla.eps", plot = volcano.p.6, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")

# sh4_renilla
sh4_renilla<-sh4_renilla%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p.7 <- ggplot(data = sh4_renilla, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  labs(title="sh4_renilla")+
  geom_text_repel(data = subset(sh4_renilla, adj.P.Val< 0.001 |adj.P.Val < 0.05 & logFC>=1 | adj.P.Val<0.05 & logFC<=-1), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)) 
volcano.p.7

ggsave("sh4_renilla.eps", plot = volcano.p.7, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")

# sh7_renilla
sh7_renilla<-sh7_renilla%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p.8 <- ggplot(data = sh7_renilla, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  labs(title="sh7_renilla")+
  geom_text_repel(data = subset(sh7_renilla, adj.P.Val< 0.000000001 |adj.P.Val < 0.05 & logFC>=2 | adj.P.Val<0.05 & logFC<=-2), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)) 
volcano.p.8 
ggsave("sh7_renilla.eps", plot = volcano.p.8, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")


# renillaAraC_renilla

renillaAraC_renilla<-renillaAraC_renilla%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p.9 <- ggplot(data = renillaAraC_renilla, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  labs(title="renillaAraC_renilla")+
  geom_text_repel(data = subset(renillaAraC_renilla, adj.P.Val< 0.0000000001 |adj.P.Val < 0.05 & logFC>=3 | adj.P.Val<0.05 & logFC<=-3), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)) 
volcano.p.9 
ggsave("renillaAraC_renilla.eps", plot = volcano.p.9, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")

# k3.1AraC_renillaAraC
k3.1AraC_renillaAraC<-k3.1AraC_renillaAraC%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p.10 <- ggplot(data = k3.1AraC_renillaAraC, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  labs(title="k3.1AraC_renillaAraC")+
  geom_text_repel(data = subset(k3.1AraC_renillaAraC, adj.P.Val< 0.001 |adj.P.Val < 0.05 & logFC>=1 | adj.P.Val<0.05 & logFC<=-1), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)) 
volcano.p.10 
ggsave("k3.1AraC_renillaAraC.eps", plot = volcano.p.10, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")

# k7AraC_renillaAraC
k7AraC_renillaAraC<-k7AraC_renillaAraC%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p.11 <- ggplot(data = k7AraC_renillaAraC, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  labs(title="k7AraC_renillaAraC")+
  geom_text_repel(data = subset(k7AraC_renillaAraC, adj.P.Val< 0.0000001 |adj.P.Val < 0.05 & logFC>=2.5 | adj.P.Val<0.05 & logFC<=-2.5), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)) 
volcano.p.11
ggsave("k7AraC_renillaAraC.eps", plot = volcano.p.11, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/plots")
############################################################################################################################

#### write tables for significant ones ####################################################################################
write.csv(sh3.1_GFP_p, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/sh3.1_GFP_p.csv")
write.csv(sh4_GFP_p, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/sh4_GFP_p.csv")
write.csv(sh7_GFP_p, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/sh7_GFP_p.csv")
write.csv(k3.1AraC_k3.1_p, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/k3.1AraC_k3.1_p.csv")
write.csv(k7AraC_k7_p, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/k7AraC_k7_p.csv")
write.csv(sh3.1_renilla_p, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/sh3.1_renilla_p.csv")
write.csv(sh4_renilla_p, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/sh4_renilla_p.csv")
write.csv(sh7_renilla_p, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/sh7_renilla_p.csv")
write.csv(renilla_GFP_p, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/renilla_GFP_p.csv")
write.csv(renillaAraC_renilla_p, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/renillaAraC_renilla_p.csv")
write.csv(k3.1AraC_renillaAraC_p, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/k3.1AraC_renillaAraC_p.csv")
write.csv(k7AraC_renillaAraC_p, file="/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/k7AraC_renillaAraC_p.csv")

###########################################################################################################################

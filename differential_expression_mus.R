########## Differential gene expression: Limma ########################################
# differential expression between 2A3 and 2A3A9
# differential expression between 4B9 and 4B9G11

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis")
load("counts.RData")
targets_mus<-read.table("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/anno/anno_mus.txt", sep="\t")

library(limma)
library(edgeR)
library(dplyr)
library(biomaRt)
library(AnnotationDbi,lib.loc="~/R/x86_64-redhat-linux-gnu-library/3.4")
source("functions.R")
#######################################################################################

##### 2A3 differential gene expression ###############################################
counts_mus_2<-counts_mus[,grepl("2A3",colnames(counts_mus))]
targets_mus_2<-targets_mus[grepl("2A3",targets_mus$sample_name),]
targets_mus_2$state2<-as.factor(as.character(targets_mus_2$state2))

# make DGE object
dge <- DGEList(counts = counts_mus_2,lib.size = colSums(counts_mus_2), samples=targets_mus_2, remove.zeros = TRUE)
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

# model matrix and contrast
design<-stats::model.matrix(~0+state,data=targets_mus_2)
colnames(design)<-c("KO","WT")

cont.matrix <- makeContrasts(
  KO_WT= KO-WT,
  levels = design)

# empirical bayes model fitting
v <- voom(dge,design = design, plot = TRUE)
fit <- lmFit(v,design)
fit2<-contrasts.fit(fit,contrasts = cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.1, trend = F)

# contrast: KO vs WT
results_2<-topTable(fit2, number = Inf, coef=1, adjust.method = "BH",confint = TRUE)
results_2_p<-results_2[results_2$adj.P.Val<=0.05,]

# get gene names
results_2<-getGeneID_mus(results_2)
results_2_p<-getGeneID_mus(results_2_p)

# decideTests
dec_tests <- decideTests(object = fit2, method = "separate", adjust.method = "BH", p.value = 0.05)
summary(dec_tests)

# volcano plot: 
results_2<-results_2%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p <- ggplot(data = results_2, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  geom_text_repel(data = subset(results_2, adj.P.Val< 0.001 |adj.P.Val < 0.05 & logFC>=1 | adj.P.Val<0.05 & logFC<=-1), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)) 
volcano.p 

ggsave("volcano_diff_exp_2A3_mus.eps", plot = volcano.p, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")

write.csv(results_2_p, file="/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/DE_2A3_KO_WT.csv")
###########################################################################################

##### 4B9 differential gene expression ###############################################
counts_mus_4<-counts_mus[,grepl("4B9",colnames(counts_mus))]
targets_mus_4<-targets_mus[grepl("4B9",targets_mus$sample_name),]
targets_mus_4$state2<-as.factor(as.character(targets_mus_4$state2))

# make DGE object
dge <- DGEList(counts = counts_mus_4,lib.size = colSums(counts_mus_4), samples=targets_mus_4, remove.zeros = TRUE)
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

# model matrix and contrast
design<-stats::model.matrix(~0+state,data=targets_mus_4)
colnames(design)<-c("KO","WT")

cont.matrix <- makeContrasts(
  KO_WT= KO-WT,
  levels = design)

# empirical bayes model fitting
v <- voom(dge,design = design, plot = TRUE)
fit <- lmFit(v,design)
fit2<-contrasts.fit(fit,contrasts = cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.1, trend = F)

# contrast: KO vs WT
results_4<-topTable(fit2, number = Inf, coef=1, adjust.method = "BH",confint = TRUE)
results_4_p<-results_4[results_4$adj.P.Val<=0.05,]

# get gene names
results_4<-getGeneID_mus(results_4)
results_4_p<-getGeneID_mus(results_4_p)

# decideTests
dec_tests <- decideTests(object = fit2, method = "separate", adjust.method = "BH", p.value = 0.05)
summary(dec_tests)

# volcano plot: 
results_4<-results_4%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p <- ggplot(data = results_4, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  geom_text_repel(data = subset(results_4, adj.P.Val <= 0.05), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)) 
volcano.p 

ggsave("volcano_diff_exp_4B9_mus.eps", plot = volcano.p, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")

write.csv(results_4_p, file="/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/DE_4B9_KO_WT.csv")


#### venn diagramm 2A3 and 4B9 ###########################################################################
overlap<-match(results_2_p$external_gene_name, results_4_p$external_gene_name)
results_overlap<-results_4_p[overlap,]
results_overlap<-na.omit(results_overlap)
# 55 genes and 9 genes the overlap is 5 genes (results_overlap)

pdf("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots/venn_2A3_4B9.pdf")
grid.newpage();
venn.plot <- draw.pairwise.venn(55, 9, 5, c("2A3", "4B9"),fill = c("blue","green"),
                                col = c("royalblue","seagreen2"),alpha=0.06, scaled=FALSE);
dev.off()

write.csv(results_overlap, file="/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/DE_overlap_2A3_4B9_mus.csv")
###########################################################################################################








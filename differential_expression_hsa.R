########## Differential gene expression: Limma ########################################
# differential expression between scramble siRNA (scr) and the two siRNAs (32, 34) for K562 and Hek293T
# differential expression between MM1 and MM6

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis")
load("counts.RData")
targets_hsa<-read.table("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/anno/anno_hsa.txt", sep="\t")

library(limma)
library(edgeR)
library(dplyr)
library(biomaRt)
library(org.Hs.eg.db) 
library(AnnotationDbi)
library(ggrepel)
library(extrafont)
source("functions.R")
loadfonts(device="postscript")
#######################################################################################

##### K562 differential gene expression ###############################################
counts_hsa_k<-counts_hsa[,grepl("K562",colnames(counts_hsa))]
targets_hsa_k<-targets_hsa[grepl("K562",targets_hsa$sample_name),]
targets_hsa_k$state2<-as.factor(as.character(targets_hsa_k$state2))

# make DGE object
dge <- DGEList(counts = counts_hsa_k,lib.size = colSums(counts_hsa_k), samples=targets_hsa_k, remove.zeros = TRUE)
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
design<-stats::model.matrix(~state2,data=targets_hsa_k)
colnames(design)<-c("Intercept","siRNA_32","siRNA_34")

cont.matrix <- makeContrasts(
  siRNA32_scr= siRNA_32,
  siRNA34_scr= siRNA_34,
  levels = design)

# empirical bayes model fitting
v <- voom(dge,design = design, plot = TRUE)
fit <- lmFit(v,design)
fit2<-contrasts.fit(fit,contrasts = cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.1, trend = F)
fit2

# contrast 1: scr_siRNA32
results_1<-topTable(fit2, number = Inf, coef=1, adjust.method = "BH",confint = TRUE)
results_1_p<-results_1[results_1$adj.P.Val<=0.05,]
# contrast 2: scr_siRNA34
results_2<-topTable(fit2, number = Inf, coef=2, adjust.method = "BH",confint = TRUE)
results_2_p<-results_2[results_2$adj.P.Val<=0.05,]

# get gene names
results_1<-getGeneID_hsa(results_1)
results_2<-getGeneID_hsa(results_2)
results_1_p<-getGeneID_hsa(results_1_p)
results_2_p<-getGeneID_hsa(results_2_p)

# decideTests
dec_tests <- decideTests(object = fit2, method = "separate", adjust.method = "BH", p.value = 0.05)
summary(dec_tests)

# venn diagram
pdf("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots/venn_K562.pdf")
limma::vennDiagram(dec_tests,include = "both",circle.col = c("royalblue","seagreen2","plum","slateblue"),cex = c(1,0.8,0.8))
dev.off()

# volcano plot: siRNA 32 vs scr
results_1<-results_1%>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p <- ggplot(data = results_1, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " adjusted p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  geom_text_repel(data = subset(results_1, adj.P.Val< 0.001 |adj.P.Val < 0.05 & logFC>=1 | adj.P.Val<0.05 & logFC<=-1 | external_gene_name =="SLC29A1"), 
                  aes(label = external_gene_name), 
                  size = 4, family= "Arial",
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 16, colour = "black",family = "Arial"),
        axis.text = element_text(size = 12, colour = "black",family = "Arial"), strip.text = element_text(size = 10, family = "Arial")) 
volcano.p 

# volcano plot: siRNA 34 vs scr
results_2<-results_2%>%mutate(threshold = ifelse(adj.P.Val <=0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p.2 <- ggplot(data = results_2, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  geom_text_repel(data = subset(results_2, adj.P.Val < 0.05 ), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)) 
volcano.p.2

ggsave("volcano_diff_exp_K562_siRNA32_scr_hsa.pdf", plot = volcano.p, width= 150, height= 150, units= "mm",  device=pdf, path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")
ggsave("volcano_diff_exp_K562_siRNA34_scr_hsa.eps", plot = volcano.p.2, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")

write.csv(results_1_p, file="/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/DE_K562_siRNA32_scr.csv")
write.csv(results_2_p, file="/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/DE_K562_siRNA34_scr.csv")
###########################################################################################

##### Hek293T differential gene expression ###############################################
counts_hsa_h<-counts_hsa[,grepl("Hek",colnames(counts_hsa))]
targets_hsa_h<-targets_hsa[grepl("Hek",targets_hsa$sample_name),]
targets_hsa_h$state2<-as.factor(as.character(targets_hsa_h$state2))

# make DGE object
dge <- DGEList(counts = counts_hsa_h,lib.size = colSums(counts_hsa_h), samples=targets_hsa_h, remove.zeros = TRUE)
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
design<-stats::model.matrix(~state2,data=targets_hsa_h)
colnames(design)<-c("Intercept","siRNA_32","siRNA_34")

cont.matrix <- makeContrasts(
  siRNA32_scr= siRNA_32,
  siRNA34_scr= siRNA_34,
  levels = design)

# empirical bayes model fitting
v <- voom(dge,design = design, plot = TRUE)
fit <- lmFit(v,design)
fit2<-contrasts.fit(fit,contrasts = cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.1, trend = F)

# contrast 1: siRNA32_scr
results_1<-topTable(fit2, number = Inf, coef=1, adjust.method = "BH",confint = TRUE)
results_1_p<-results_1[results_1$adj.P.Val<=0.05,]
# contrast 2: siRNA34_scr
results_2<-topTable(fit2, number = Inf, coef=2, adjust.method = "BH",confint = TRUE)
results_2_p<-results_2[results_2$adj.P.Val<=0.05,]

# get gene names
results_1<-getGeneID_hsa(results_1)
results_2<-getGeneID_hsa(results_2)
results_1_p<-getGeneID_hsa(results_1_p)
results_2_p<-getGeneID_hsa(results_2_p)

# decideTests
dec_tests <- decideTests(object = fit2, method = "separate", adjust.method = "BH", p.value = 0.05)
summary(dec_tests)

# venn diagram
pdf("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots/venn_Hek293T.pdf")
limma::vennDiagram(dec_tests,include = "both",circle.col = c("royalblue","seagreen2","plum","slateblue"),cex = c(1,0.8,0.8))
dev.off()

# volcano plot: siRNA34 vs scr
results_2<-results_2%>%mutate(threshold = ifelse(adj.P.Val <=0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C")))

volcano.p <- ggplot(data = results_2, aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  geom_text_repel(data = subset(results_2, adj.P.Val < 0.05), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)) 
volcano.p

ggsave("volcano_diff_exp_HEK_siRNA32_scr_hsa.eps", plot = volcano.p, width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/plots")

write.csv(results_1_p, file="/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/DE_Hek293T_siRNA32_scr.csv")
write.csv(results_2_p, file="/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/DE_Hek293T_siRNA34_scr.csv")
###########################################################################################

##### MM1 vs MM6 differential gene expression ###############################################
counts_hsa_m<-counts_hsa[,grepl("MM",colnames(counts_hsa))]
targets_hsa_m<-targets_hsa[grepl("MM",targets_hsa$sample_name),]
targets_hsa_m$cell_line<-as.factor(as.character(targets_hsa_m$cell_line))

# make DGE object
dge <- DGEList(counts = counts_hsa_m,lib.size = colSums(counts_hsa_m), samples=targets_hsa_m, remove.zeros = TRUE)
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
design<-stats::model.matrix(~0+cell_line,data=targets_hsa_m)
colnames(design)<-c("MM1","MM6")

cont.matrix <- makeContrasts(
  MM1_MM6= MM1-MM6,
  levels = design)

# empirical bayes model fitting
v <- voom(dge,design = design, plot = TRUE)
fit <- lmFit(v,design)
fit2<-contrasts.fit(fit,contrasts = cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.1, trend = F)

# contrast: MM1 vs MM6
results<-topTable(fit2, number = Inf, coef=1, adjust.method = "BH",confint = TRUE)
results_p<-results[results$adj.P.Val<=0.05,]

# get gene names
results<-getGeneID_hsa(results)
results_p<-getGeneID_hsa(results)


# decideTests
dec_tests <- decideTests(object = fit2, method = "separate", adjust.method = "BH", p.value = 0.05)
summary(dec_tests)

# subset results with logFC higher than 2 or smaller than -2
results_logFC2<-results[results$logFC>=2|results$logFC<=-2,]

# subset results with logFC higher than 4 or smaller than -4
results_logFC4<-results[results$logFC>=4|results$logFC<=-4,]

write.csv(results, file="/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/DE_MM1_MM6.csv")
write.csv(results_logFC2, file="/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/DE_MM1_MM6_logFC2.csv")
write.csv(results_logFC4, file="/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/DE_MM1_MM6_logFC4.csv")

# differential expression of NT5C2
NT5C2<-results[grepl("NT5C2", results$external_gene_name),]
write.csv(NT5C2, file="/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/DE_MM1_MM6_NT5C2.csv")
###########################################################################################

# differential expression of SPARC
SPARC<-results[grepl("SPARC", results$external_gene_name),]
write.csv(SPARC, file="/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis/tables/DE_MM1_MM6_SPARC.csv")
###########################################################################################






########## Differential gene expression: Limma ########################################

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018")
load("counts_k562_march.RData")

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


#####################################################################################################################

#### model matrix and contrast: contrasts of sh RNAs ###############################################################
design<-stats::model.matrix(~0+sample,data=anno)
colnames(design)<-gsub("sample","",colnames(design))

cont.matrix <- makeContrasts(
  A4_nativ = A4_No - nativ_No,
  F7_nativ = F7_No - nativ_No,
  C5_nativ = C5_No - nativ_No,
  A4_WT2 = A4_No - WT2_No,
  F7_WT2 = F7_No - WT2_No,
  C5_WT2 = C5_No - WT2_No,
  F7.AraC24h_nativ.AraC24h = F7_AraC24h - nativ_AraC24h,
  F7.AraC48h_nativ.AraC48h = F7_AraC48h - nativ_AraC48h,
  F7.AraC72h_nativ.AraC72h = F7_AraC72h - nativ_AraC72h,
  nativ_AraC24h_nativ = nativ_AraC24h - nativ_No,
  nativ_AraC48h_nativ = nativ_AraC48h - nativ_No,
  nativ_AraC72h_nativ = nativ_AraC72h - nativ_No,
  F7.AraC24h_F7 = F7_AraC24h - F7_No,
  F7.AraC48h_F7 = F7_AraC48h - F7_No,
  F7.AraC72h_F7 = F7_AraC72h - F7_No,
  levels = design)

# empirical bayes model fitting
v <- voom(dge,design = design, plot = TRUE)
fit <- lmFit(v,design)
fit2<-contrasts.fit(fit,contrasts = cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.1, trend = F)

# get the contrasts
# A4_nativ<-topTable(fit2, number = Inf, coef=j, adjust.method = "BH",confint = TRUE)
# A4_nativ_p<-A4_nativ[A4_nativ$adj.P.Val<=0.05,]
# get gene names with getGeneID_hsa function
# A4_nativ<-getGeneID_hsa(A4_nativ)
# A4_nativ_p<-getGeneID_hsa(A4_nativ_p)

j<-0
for (i in colnames(cont.matrix)){
j=j+1
  assign(paste(i,"_g",sep = ""), getGeneID_hsa(assign(paste(i), topTable(fit2, number = Inf, coef=j, adjust.method = "BH",confint = TRUE))))
  assign(paste(i,"_p",sep = ""),getGeneID_hsa(assign(paste(i), topTable(fit2, number = Inf, coef=j, adjust.method = "BH",confint = TRUE))[assign(paste(i), topTable(fit2, number = Inf, coef=j, adjust.method = "BH",confint = TRUE))$adj.P.Val<=0.05,]))
}

#######################################################################################################################

save(list=ls(),
     file="diff.RData");
load("diff.RData")
#### decideTests ######################################################################################################
dec_tests <- decideTests(object = fit2, method = "separate", adjust.method = "BH", p.value = 0.05)
summary(dec_tests)
#######################################################################################################################

#### Venn diagrams ####################################################################################################
nativ<-dec_tests[,10:12]
clones_nativ<-dec_tests[,1:3]
clones_WT2<-dec_tests[,4:6]
F7<-dec_tests[,7:9]

tests<-list(nativ=nativ,clones_nativ=clones_nativ,clones_WT2=clones_WT2,F7=F7)

pdf("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/plots/venn.pdf")
for (i in tests){
limma::vennDiagram(i, include = "both",circle.col = c("royalblue","seagreen2","plum","slateblue"),cex = c(1,0.8,0.8))
}
dev.off()
####################################################################################################################

#### which genes do overlap ######################################################################################################

overlap <- VennDiagram::calculate.overlap(list(rownames(A4_WT2_g[A4_WT2_g$adj.P.Val<=0.05,]),
                                               rownames(F7_WT2_g[F7_WT2_g$adj.P.Val<=0.05,]),
                                               rownames(C5_WT2_g[C5_WT2_g$adj.P.Val<=0.05,])))

# a1-a7: circles in venn diagram (start counting from the left to right, then second row, left to right etc.)
a5<-A4_WT2_g %>% dplyr::filter(rownames(A4_WT2_g)%in% overlap$a5)
a2<-A4_WT2_g %>% dplyr::filter(rownames(A4_WT2_g) %in% overlap$a2)

write.csv(a5, file="/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/overlap_A4WT2_F7WT2_C5WT2.csv")
write.csv(a2, file="/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/overlap_A4WT2_F7WT2.csv")

#### 
decidetests global method



## looop for getting the overlap from all comparisions 
k<-0
for (i in tests){
  for (j in colnames(i)){
    k=k+1
    print(paste(k))
    assign(paste("row",k,sep = "_"), rownames(get(j)))
    list<-list()
  }
}



#######################################################

##### volcano plot #######################################################################################################
for (i in colnames(cont.matrix)){
  assign(paste(i,"_g",sep = ""),get(paste(i,"_g",sep = "")) %>%mutate(threshold = ifelse(adj.P.Val <= 0.05,"A" , ifelse(adj.P.Val>0.05, "B", "C"))))

assign(paste("volcano_", i, sep = ""), ggplot(data = get(paste(i,"_g",sep = "")), aes(x = logFC, y = -log10(adj.P.Val))) + 
  geom_point(aes(colour = threshold)) +
  ylab(expression(paste(-log[10], " p-value"))) +
  xlab(expression(paste(log[2],"FC")))+
  labs(title=paste(i))+
  geom_text_repel(data = subset(get(paste(i,"_g",sep = "")), adj.P.Val < 0.05 & logFC>=2 | adj.P.Val<0.05 & logFC<=-2), 
                  aes(label = external_gene_name), 
                  size = 3, 
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"))+
  theme_bw() + geom_hline(yintercept = -log10(0.05), linetype = "dashed",colour = "grey") + 
  theme(legend.position = "none", axis.title = element_text(size = 10, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"), strip.text = element_text(size = 6)))
paste("volcano_", i, ".eps",sep = "")
#ggsave(paste("volcano_", i, ".eps",sep = ""), plot = get(paste("volcano_", i, sep = "")), width= 150, height= 150, units= "mm", device = "eps",path = "/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/plots")
}

#### write tables for diff exp genes (significant) ####################################################################################
for (i in colnames(cont.matrix)){
write.csv(get(paste(i,"_p",sep = "")), file=paste("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/",i,"_p.csv",sep = ""))
}
###########################################################################################################################

#### write whole table of genes ####################################################################################
for (i in colnames(cont.matrix)){
  write.csv(get(paste(i,"_g",sep = "")), file=paste("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/DE_complete/",i,"_g.csv",sep = ""))
}
###########################################################################################################################


###### Ko F7 with AraC against nativ with ARAC ###########################################################################
cont.matrix <- makeContrasts(
  F7.AraC24h_nativ.AraC24h = F7_AraC24h - nativ_AraC24h,
  F7.AraC48h_nativ.AraC48h = F7_AraC48h - nativ_AraC48h,
  F7.AraC72h_nativ.AraC72h = F7_AraC72h - nativ_AraC72h,
  levels = design)

v <- voom(dge,design = design, plot = TRUE)
fit <- lmFit(v,design)
fit2<-contrasts.fit(fit,contrasts = cont.matrix)
fit2 <- eBayes(fit2, proportion = 0.1, trend = F)

result<-topTableF(fit2, number = Inf, adjust.method = "BH")
# adjp.value for all samples
# logFC per contrast

dec_tests <- decideTests(object = fit2, method = "global", adjust.method = "BH", p.value = 0.05)
summary(dec_tests)

data<-as.data.frame(dec_tests@.Data)
up<-data[data$F7.AraC24h_nativ.AraC24h==1 & data$F7.AraC48h_nativ.AraC48h==1 & data$F7.AraC72h_nativ.AraC72h==1,]
down<-data[data$F7.AraC24h_nativ.AraC24h==-1 & data$F7.AraC48h_nativ.AraC48h==-1 & data$F7.AraC72h_nativ.AraC72h==-1,]

up_gene<-rownames(up)
down_gene<-rownames(down)

library(biomaRt)
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
listAttributes(ensembl)
listFilters(ensembl)

go_up = getBM(attributes = c("external_gene_name","ensembl_gene_id"), 
           filters = c("ensembl_gene_id"), 
           values = c(up_gene), 
           mart = ensembl)

go_down = getBM(attributes = c("external_gene_name","ensembl_gene_id"), 
              filters = c("ensembl_gene_id"), 
              values = c(down_gene), 
              mart = ensembl)


write.csv(go_up, file = "/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/F7Arac_nativArac_up.csv")
write.csv(go_down, file = "/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/F7Arac_nativArac_down.csv")

#### GSEA Gene Set Enrichment Analysis gseGO
library(clusterProfiler)
library(org.Hs.eg.db)

#ranked list of genes
genes<-result$adj.P.Val
names(genes)<-rownames(result)
genes<-sort(genes, decreasing = TRUE)

# gene set enrichment
gsea <- gseGO(gene=genes,OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP",
              pAdjustMethod = "BH", pvalueCutoff  = 0.1, by="fgsea")
gsea_result<-gsea@result

# no terms enriched in GSEA!!

### GO Analysis

genes<-c(go_down$ensembl_gene_id, go_up$ensembl_gene_id)

keytypes(org.Hs.eg.db)


ego <- enrichGO(gene=genes,OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,readable = TRUE)
ego_result<-ego@result

# simplify result, combine closely related terms
ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
ego2_result<-ego2@result

# significant results
ego2_result_sign<-ego2_result[ego2_result$p.adjust<0.05,]


write.csv(ego2_result_sign, file ="/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/GO_Enrichment.csv") 



d<-dotplot(ego2, showCategory=20, font.size = 10)
d
ggsave("GO_Enrichment_dotplot.pdf", plot= d, device = "pdf", path = "/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/plots",
       units="cm", width=20, height=20)


e<-emapplot(ego, showCategory = 20)
e
ggsave("GO_Enrichment_emapplot.pdf", plot= e, device = "pdf", path = "/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/plots",
       units="cm", width=20, height=20)


# vector of up or down for cnetplot
# mean of lfc of all contrasts
result_sub<-result[result$adj.P.Val<0.05,c(1:3)]
lfc<-result$adj.P.Val

lfc<-apply(result,1,mean)


length(lfc)

cnetplot(ego, showCategory = 5, foldChange=2^lfc,vertex.label.cex = 1.2 )+
  labs(title = "Temperatures\n", fill="bla") 





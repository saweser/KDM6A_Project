rm(list=ls());  # empty workspace
setwd("/Users/sabrina/Documents/")

library(clusterProfiler)
library(org.Hs.eg.db)
genes<-read.csv("/Users/sabrina/Downloads/A4_WT2_p.csv")

keytypes(org.Hs.eg.db)

ego <- enrichGO(gene=genes$X,OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,readable = TRUE)

barplot(ego, drop=TRUE, showCategory=12)
# showCategory = how many bars


View(ego[1:100,])

a<-ego@result
barplot(as.matrix(a), showCategory=12)

goplot(ego)


# DOSE
genes_entrez = bitr(genes$X, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

edo <- enrichDO(gene=genes_entrez$ENTREZID, ont = "DO",
                pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,readable = FALSE)
b<-edo@result

##### GSEA #####################################################################################################################
rm(list=ls());  # empty workspace
setwd("/Users/sabrina/Documents/")

library(clusterProfiler)
library(org.Hs.eg.db)


# make a list of all files to be processed
file.list=as.list(list.files("/Users/sabrina/Downloads/DE_complete", pattern=".csv",full.names=TRUE))
files=lapply(file.list, read.csv)
names(files)<-gsub("/Users/sabrina/Downloads/DE_complete/","",file.list)

# add rank, logFC*1-adj.P.Val
for(i in 1:length(files)){
  files[[i]]$rank <- files[[i]]$logFC*(1-(files[[i]]$adj.P.Val))
}

# chose file to do GSEA, blah= number of list object
for (i in 1:15){
  blah<-i
  genes_gsea<-as.numeric(files[[blah]]$rank)
  names(genes_gsea)<-as.character(files[[blah]]$X)
  
  # sort genes according to rank
  genes_gsea<-sort(genes_gsea, decreasing = TRUE)
  
  # gene set enrichment
  gsea <- gseGO(gene=genes_gsea,OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff  = 0.05, by="fgsea")
  gsea_result<-gsea@result
  
  # save list of GOs
  name<-gsub("_g.csv","",names(files[blah]))
  write.csv(gsea_result, file = paste("/Users/sabrina/Documents/plots/", name,"_gsea.csv",sep = ""))
}

# plot Gene Ontology (GO) of interest
gseaplot(gsea, geneSetID = "GO:0018345")
################################################################################################################################







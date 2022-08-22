rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017")
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)

anno<-read_excel("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/anno_k562.xlsx")
###############################################################################################################################
anno<-as.data.frame(anno)
rownames(anno)<-anno$sample_name

# add column treatment
anno$treatment<-rownames(anno)
anno$treatment<-gsub(pattern = "K562_", replacement = "",anno$treatment)
anno$treatment<-gsub(pattern = "_6617",replacement = "",anno$treatment)
anno$treatment<-gsub(pattern = "_7617",replacement = "",anno$treatment)
anno$treatment<-gsub(pattern = "_9617",replacement = "",anno$treatment)
anno$treatment<-gsub(pattern = "_11617",replacement = "",anno$treatment)
anno$treatment<-gsub(pattern = "_12617",replacement = "",anno$treatment)
anno$treatment<-gsub(pattern = "_15617",replacement = "",anno$treatment)
anno$treatment<-gsub(pattern = "_17617",replacement = "",anno$treatment)
anno$treatment<-gsub(pattern = "Nr.1",replacement = "",anno$treatment)
anno$treatment<-gsub(pattern = "Nr.2",replacement = "",anno$treatment)

# add column shRNA_2
anno$shRNA_2<-anno$treatment
anno$shRNA_2<-gsub(pattern = "_AraC",replacement = "",anno$shRNA_2)

# add column shRNA
anno$shRNA<-anno$shRNA_2
anno$shRNA<-gsub(pattern = "_.*",replacement = "",anno$shRNA)

# add AraC
anno$AraC<-anno$treatment
anno$AraC<-gsub(pattern = ".*_",replacement = "",anno$AraC)
anno$AraC<-gsub(pattern = "shRenilla",replacement = "No",anno$AraC)
anno$AraC<-gsub(pattern = "shGFP",replacement = "No",anno$AraC)
anno$AraC<-gsub(pattern = "[[:digit:]]",replacement = "No",anno$AraC)
anno$AraC<-gsub(pattern = "No.No",replacement = "No",anno$AraC)
anno$AraC<-gsub(pattern = "AraC",replacement = "Yes",anno$AraC)

# add column date
anno$date<-gsub(pattern = ".*_",replacement = "",anno$sample_name)

# add cell line
anno$cell_line<-gsub(pattern = "_.*",replacement = "",anno$sample_name)

# make sure that the rownames of annotation and the colanmes of counts are in the same oder
# check with colnamescounts==rownames(anno)
anno<-anno[order(anno$sample_name),]


# save the annotation table
write.table(anno, file = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/anno.txt", col.names = TRUE, row.names = TRUE)




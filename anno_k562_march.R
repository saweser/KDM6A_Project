rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018")
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)

anno<-read.table("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/anno_k562_march.txt", header = T)
counts<-read.table("/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/counts_k562_march2018.txt",header = T)
###############################################################################################################################
rownames(anno)<-anno$name

# add column state
anno$state<-rownames(anno)
anno$state<-gsub(pattern = "K562_", replacement = "",anno$state)
anno$state<-gsub(pattern = "[+]", replacement = "_",anno$state)
anno$state<-gsub(pattern = "_.*", replacement = "",anno$state)

# add column treatment
anno$treatment<-rownames(anno)
anno$treatment<-gsub(pattern = ".*[+]", replacement = "",anno$treatment)
anno$treatment<-gsub(pattern = "_.*", replacement = "",anno$treatment)
anno$treatment<-gsub(pattern = "K562", replacement = "No",anno$treatment)

# add column sample
anno<-anno %>% unite(sample, c("state","treatment"), sep = "_", remove = FALSE)

# add column date
anno$date<-gsub(pattern = ".*_",replacement = "",anno$name)

# save the annotation table
write.table(anno, file = "/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/tables/anno.txt", col.names = TRUE, row.names = TRUE)




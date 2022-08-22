##### count tables ##############################################################################################################
# make the count tables for human and mouse

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis")
AllCounts<-readRDS("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/zUMIs_output/expression/Sophie.dgecounts.rds")
sample.names<-read.table("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/anno/SCRBseq_96_setA_Sample_names.txt",sep = "\t")
#################################################################################################################################

##### downsampled UMI counts
names(AllCounts$exons$downsampled$downsampled_8374190$umicounts_downsampled)
counts<-AllCounts$exons$downsampled$downsampled_8374190$umicounts_downsampled

#### change column names
sample.names<-sample.names[,2:4]
sample.names<-na.omit(sample.names)
colnames(sample.names)<-c("well","barcode","sample.name")
sample.names$barcode<-as.character(sample.names$barcode)

identical(colnames(counts),sample.names$barcode)
all.equal(colnames(counts),sample.names$barcode)

########
counts$gene<-rownames(counts)
counts.long <- reshape2::melt(counts)
colnames(counts.long)<-c("gene","barcode","value")
counts.long<-left_join(counts.long,sample.names,by="barcode")
counts.long$barcode<-NULL
counts.long$well<-NULL
counts<-dcast(counts.long,gene~sample.name,value.var="value")
rownames(counts)<-counts$gene
counts$gene<-NULL
#################

#### subset human and mouse
# add column for human and mouse
sample.names$genome<-""
for (i in 1:nrow(sample.names)){                            
  if (grepl("2A3",sample.names[i, "sample.name"])){           
    sample.names[i, "genome"]<-"mouse"                   
  } else if(grepl("4B9",sample.names[i, "sample.name"])){
    sample.names[i, "genome"]<-"mouse"
  }else{
    sample.names[i, "genome"]<-"human"
  }
  }  

# mouse genes are labelled as ENSMUSG!
# human
human<-sample.names[sample.names$genome=="human",]
human<-as.character(human$sample.name)
counts_hsa<-counts[,human]
counts_hsa<-counts_hsa[grepl("ENSG",rownames(counts_hsa)),]

# mouse
mouse<-sample.names[sample.names$genome=="mouse",]
mouse<-as.character(mouse$sample.name)
counts_mus<-counts[,mouse]
counts_mus<-counts_mus[grepl("ENSMUSG",rownames(counts_mus)),]


##### save objects for downstream analysis ###############################################################
save(list=c("counts_hsa","counts_mus"),
     file="counts.RData");
##############################################################################################################



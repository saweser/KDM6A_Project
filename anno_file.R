########## make the annotation file #######################################################################

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/analysis")
anno<-read.table("/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/anno/SCRBseq_96_setA_Sample_names.txt",sep = "\t")
###########################################################################################################


anno<-na.omit(anno)
colnames(anno)<-c("set","well","barcode","sample_name")
rownames(anno)<-anno$sample_name

# add cell line
anno$cell_line<-""
for (i in 1:nrow(anno)){                            
  if (grepl("MM1",anno[i, "sample_name"])){           
    anno[i, "cell_line"]<-"MM1"                   
  } else if(grepl("MM6",anno[i, "sample_name"])){
    anno[i, "cell_line"]<-"MM6"
  }else if(grepl("Hek",anno[i, "sample_name"])){
    anno[i, "cell_line"]<-"HEK"
  }else if(grepl("4B9",anno[i, "sample_name"])){
    anno[i, "cell_line"]<-"4B9"
  }else if(grepl("K562",anno[i, "sample_name"])){
    anno[i, "cell_line"]<-"K562"
  }else if(grepl("2A3",anno[i, "sample_name"])){
    anno[i, "cell_line"]<-"2A3"
  }
}  

# add mutation state
anno$state<-""
for (i in 1:nrow(anno)){                            
  if (grepl("MM1",anno[i, "sample_name"])){           
    anno[i, "state"]<-"WT"  
  }else if(grepl("2A3A9",anno[i, "sample_name"])){
    anno[i, "state"]<-"KO"
  }else if(grepl("4B9G11",anno[i, "sample_name"])){
    anno[i, "state"]<-"KO"
  }else if(grepl("2A3A9",anno[i, "sample_name"])){
    anno[i, "state"]<-"KO"
  }else if(grepl("Heksi",anno[i, "sample_name"])){
    anno[i, "state"]<-"siRNA"
  }else if(grepl("Hekscr",anno[i, "sample_name"])){
    anno[i, "state"]<-"scr"
  }else if(grepl("K562scr",anno[i, "sample_name"])){
    anno[i, "state"]<-"scr"
  }else if(grepl("K562si",anno[i, "sample_name"])){
    anno[i, "state"]<-"siRNA"
  }else{
    anno[i, "state"]<-"WT"
  }
}  

# add 2nd mutation state to distinguish between siRNA 32 and 34
anno$state2<-""
for (i in 1:nrow(anno)){                            
  if (grepl("MM1",anno[i, "sample_name"])){           
    anno[i, "state2"]<-"WT"  
  }else if(grepl("2A3A9",anno[i, "sample_name"])){
    anno[i, "state2"]<-"KO"
  }else if(grepl("4B9G11",anno[i, "sample_name"])){
    anno[i, "state2"]<-"KO"
  }else if(grepl("2A3A9",anno[i, "sample_name"])){
    anno[i, "state2"]<-"KO"
  }else if(grepl("Heksi32",anno[i, "sample_name"])){
    anno[i, "state2"]<-"siRNA_32"
  }else if(grepl("Heksi34",anno[i, "sample_name"])){
    anno[i, "state2"]<-"siRNA_34"
  }else if(grepl("Hekscr",anno[i, "sample_name"])){
    anno[i, "state2"]<-"scr"
  }else if(grepl("K562scr",anno[i, "sample_name"])){
    anno[i, "state2"]<-"scr"
  }else if(grepl("K562si32",anno[i, "sample_name"])){
    anno[i, "state2"]<-"siRNA_32"
  }else if(grepl("K562si34",anno[i, "sample_name"])){
    anno[i, "state2"]<-"siRNA_34"
  }else{
    anno[i, "state2"]<-"WT"
  }
}  



# add genome
anno$genome<-""
for (i in 1:nrow(anno)){                            
  if (grepl("2A3",anno[i, "sample_name"])){           
    anno[i, "genome"]<-"mouse"  
  }else if(grepl("4B9",anno[i, "sample_name"])){
    anno[i, "genome"]<-"mouse"
  }else{
    anno[i, "genome"]<-"human"
  }
}  


# add date RNA isolated
anno$name<-as.character(anno$sample_name)
anno$name[11]<-"Heksi34_026"
anno$name[44]<-"Heksi34_026"
name<-as.character(anno$name)
sbt = strsplit(name,split="_")
ncol = max(sapply(sbt,length))                        # looks at length of each row of sbt
tempsplit <- as.data.frame(lapply(1:ncol,function(i)sapply(sbt,"[",i))) 
colnames(tempsplit)<-c("name","date")
anno$date<-tempsplit$date
anno$date<-as.character(anno$date)
anno$name<-NULL




anno_hsa<-anno[anno$genome=="human",]
anno_mus<-anno[anno$genome=="mouse",]

write.table(anno, "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/anno/anno.txt", sep="\t")
write.table(anno_hsa, "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/anno/anno_hsa.txt", sep="\t")
write.table(anno_mus, "/data/htp/A07/RNA_Seq/Sophie/mESC_K562_Hek_siRNA_Jul2017/anno/anno_mus.txt", sep="\t")






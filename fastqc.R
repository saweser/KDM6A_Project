rm(list=ls());  # empty workspace
setwd("/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017")

library(fastqcr)
fastqc(fq.dir = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs", 
qc.dir = "/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/fastqc", threads = 4)

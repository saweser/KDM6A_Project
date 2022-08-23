#!/bin/bash
#SBATCH --error=fastqc.%J.err
#SBATCH --output=fastqc.%J.out
#SBATCH --workdir=/data/htp/A07/RNA_Seq/Sophie/k562_March_2018
#SBATCH --mem=100000



fastqc --threads 10 \
--outdir /data/htp/A07/RNA_Seq/Sophie/k562_March_2018/fastqc \
/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/fastq_files/180611_L183_0419_ACCGENANXX/k562_march_R1.fastq.gz \
/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/fastq_files/180611_L183_0419_ACCGENANXX/k562_march_R2.fastq.gz \
/data/htp/A07/RNA_Seq/Sophie/k562_March_2018/fastq_files/180611_L183_0419_ACCGENANXX/k562_march_R3.fastq.gz \

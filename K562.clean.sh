#!/bin/bash
#SBATCH -n 1
#SBATCH --error=clean.%J.err
#SBATCH --output=clean.%J.out
#SBATCH --workdir=/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs
#SBATCH --dependency=afterok:12546
srun rm /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.Aligned.out.bam /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.aligned.sorted.bam.in /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.aligned.sorted.bam.ex /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.barcodelist.filtered.sam
srun mv /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.barcoderead.filtered.fastq.gz /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/zUMIs_output/filtered_fastq/
srun mv /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.cdnaread.filtered.fastq.gz /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/zUMIs_output/filtered_fastq/

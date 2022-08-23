#!/bin/bash
#SBATCH -n 1
#SBATCH --error=clean.%J.err
#SBATCH --output=clean.%J.out
#SBATCH --workdir=/data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie/
#SBATCH --dependency=afterok:10257
srun rm /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.Aligned.out.bam /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.aligned.sorted.bam.in /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.aligned.sorted.bam.ex /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.barcodelist.filtered.sam
srun mv /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.barcoderead.filtered.fastq.gz /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//zUMIs_output/filtered_fastq/
srun mv /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.cdnaread.filtered.fastq.gz /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//zUMIs_output/filtered_fastq/

#!/bin/bash
#SBATCH -n 1
#SBATCH --error=unmap.%J.err
#SBATCH --output=unmap.%J.out
#SBATCH --cpus-per-task=20
#SBATCH --workdir=/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs
#SBATCH --dependency=afterok:12543
#SBATCH --mem=40
srun samtools sort -n -O sam -T tmp.K562 -@ 20 -m 2G -o /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.barcodelist.filtered.sort.sam /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.barcodelist.filtered.sam

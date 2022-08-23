#!/bin/bash
#SBATCH --error=fastqc.%J.err
#SBATCH --output=fastqc.%J.out
#SBATCH --workdir=/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=8



srun  Rscript /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/fastqc.R

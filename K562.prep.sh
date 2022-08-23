#!/bin/bash
#SBATCH -n 1
#SBATCH --error=prep.%J.err
#SBATCH --output=prep.%J.out
#SBATCH --cpus-per-task=20
#SBATCH --workdir=/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs
srun perl /data/sfb1243cs/htp/A14/zUMIs/zUMIs/fqfilter.pl /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562_Bria.R1.fq.gz /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562_Bria.R3.fq.gz 1 20 1 20 1-6 7-16 20 K562 /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs

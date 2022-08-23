#!/bin/bash
#SBATCH -n 1
#SBATCH --error=prep.%J.err
#SBATCH --output=prep.%J.out
#SBATCH --cpus-per-task=8
#SBATCH --workdir=/data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie/
srun perl /data/htp/A14/zUMIs/zUMIs/fqfilter.pl /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie/R1.fastq /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie/R3.fastq 2 20 3 20 1-6 7-16 8 Sophie /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie/

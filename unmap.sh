#!/bin/bash
#SBATCH -n 1
#SBATCH --error=unmap.%J.err
#SBATCH --output=unmap.%J.out
#SBATCH --cpus-per-task=8
#SBATCH --workdir=/data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie/
#SBATCH --dependency=afterok:10254
#SBATCH --mem=16
srun samtools sort -n -O sam -T tmp.Sophie -@ 8 -m 2G -o /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.barcodelist.filtered.sort.sam /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.barcodelist.filtered.sam

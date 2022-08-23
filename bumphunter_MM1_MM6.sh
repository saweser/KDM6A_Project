#!/bin/bash
#SBATCH --error=bumphunter.%J.err
#SBATCH --output=bumphunter.%J.out
#SBATCH --workdir=/data/htp/A07/MM1_vs_MM6
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=4



srun  Rscript /data/htp/A07/MM1_vs_MM6/bumphunter_MM1_MM6.R

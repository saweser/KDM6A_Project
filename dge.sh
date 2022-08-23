#!/bin/bash
#SBATCH -n 1
#SBATCH --error=dge.%J.err
#SBATCH --output=dge.%J.out
#SBATCH --workdir=/data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie/
#SBATCH --dependency=afterok:10255:10256
#SBATCH --mem=50000
#SBATCH --cpus-per-task=8
ln -s -f /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.aligned.sorted.bam /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.aligned.sorted.bam.in
ln -s -f /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.aligned.sorted.bam /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.aligned.sorted.bam.ex
srun --chdir=/data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie/ Rscript /data/htp/A14/zUMIs/zUMIs/zUMIs-dge.R --gtf /data/htp/A14/JohannesPDX/HsapMmus/Homo_sapiens.GRCh38.84_Mus_musculus.GRCm38.85.gtf --abam /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.aligned.sorted.bam --ubam /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.barcodelist.filtered.sort.sam --barcodefile NA --out /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie/ --sn Sophie --cores 8 --strandedness 0 --bcstart 1 --bcend 6 --umistart 7 --umiend 16 --subsamp 0
srun --chdir=/data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie/ Rscript /data/htp/A14/zUMIs/zUMIs/zUMIs-stats.R --out /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie/ --sn Sophie

#!/bin/bash
#SBATCH -n 1
#SBATCH --error=dge.%J.err
#SBATCH --output=dge.%J.out
#SBATCH --workdir=/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs
#SBATCH --dependency=afterok:12544:12545
#SBATCH --mem=50000
#SBATCH --cpus-per-task=20
ln -s -f /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.aligned.sorted.bam /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.aligned.sorted.bam.in
ln -s -f /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.aligned.sorted.bam /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.aligned.sorted.bam.ex
srun --chdir=/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs Rscript /data/sfb1243cs/htp/A14/zUMIs/zUMIs/zUMIs-dge.R --gtf /data/ngs/genomes/Human/hg38/Homo_sapiens.GRCh38.84_chrsNamesUCSC.gtf --abam /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.aligned.sorted.bam --ubam /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.barcodelist.filtered.sort.sam --barcodefile /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/tables/barcodes.txt --out /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs --sn K562 --cores 20 --strandedness 0 --bcstart 1 --bcend 6 --umistart 7 --umiend 16 --subsamp 0
srun --chdir=/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs Rscript /data/sfb1243cs/htp/A14/zUMIs/zUMIs/zUMIs-stats.R --out /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs --sn K562

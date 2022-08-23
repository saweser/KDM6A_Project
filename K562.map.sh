#!/bin/bash
#SBATCH -n 1
#SBATCH --error=map.%J.err
#SBATCH --output=map.%J.out
#SBATCH --cpus-per-task=20
#SBATCH --workdir=/data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs
#SBATCH --dependency=afterok:12543
#SBATCH --mem=28805
srun STAR --genomeDir /data/ngs/genomes/Human/hg38/STAR5idx_noGTF --runThreadN 20 --readFilesCommand zcat --sjdbGTFfile /data/ngs/genomes/Human/hg38/Homo_sapiens.GRCh38.84_chrsNamesUCSC.gtf --outFileNamePrefix /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562. --outSAMtype BAM Unsorted --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --sjdbOverhang 49 --twopassMode Basic --readFilesIn /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.cdnaread.filtered.fastq.gz 
srun samtools sort -n -O bam -T temp.K562 -@ 20 -m 2G -o /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.aligned.sorted.bam /data/htp/A07/RNA_Seq/Sophie/k562_Sept2017/run_zUMIs/K562.Aligned.out.bam

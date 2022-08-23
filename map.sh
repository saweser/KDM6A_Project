#!/bin/bash
#SBATCH -n 1
#SBATCH --error=map.%J.err
#SBATCH --output=map.%J.out
#SBATCH --cpus-per-task=8
#SBATCH --workdir=/data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie/
#SBATCH --dependency=afterok:10254
#SBATCH --mem=55811
srun STAR --genomeDir /data/htp/A14/JohannesPDX/HsapMmus/STAR5idx_noGTF/ --runThreadN 8 --readFilesCommand zcat --sjdbGTFfile /data/htp/A14/JohannesPDX/HsapMmus/Homo_sapiens.GRCh38.84_Mus_musculus.GRCm38.85.gtf --outFileNamePrefix /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie. --outSAMtype BAM Unsorted --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --sjdbOverhang 18 --twopassMode Basic --readFilesIn /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.cdnaread.filtered.fastq.gz 
srun samtools sort -n -O bam -T temp.Sophie -@ 8 -m 2G -o /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.aligned.sorted.bam /data/htp/A14/UMIpool_Jul2017/Spiekermann_Sophie//Sophie.Aligned.out.bam

#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00

nthreads=12

x=("CB14_Het" "CB42_Het" "CB43_Het" "CB20_KO" "CB38_KO" "CB41_KO")

module load STAR/2.7.3a-GCC-9.3.0
module load SAMtools/1.16.1-GCC-11.3.0

for x in "${x[@]}"; do
	STAR --genomeDir ../../Master/meta_mice/STAR_mm10/ \
		--runThreadN 12 \
		--readFilesCommand zcat \
		--readFilesIn input/${x}_1.fq.gz input/${x}_2.fq.gz \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix output/STAR/raw/${x}_
    samtools index output/STAR/raw/${x}_Aligned.sortedByCoord.out.bam
done


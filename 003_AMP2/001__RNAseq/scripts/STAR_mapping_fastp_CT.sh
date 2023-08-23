#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00

nthreads=12

x=("CT14_Het" "CT42_Het" "CT43_Het" "CT20_KO" "CT38_KO" "CT41_KO")

module load STAR/2.7.3a-GCC-9.3.0
module load SAMtools/1.16.1-GCC-11.3.0

for x in "${x[@]}"; do
	STAR --genomeDir ../../Master/meta_mice/STAR_mm10/ \
		--runThreadN 12 \
		--readFilesCommand zcat \
		--readFilesIn output/fastp/${x}_1.fq.gz output/fastp/${x}_2.fq.gz \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix output/STAR/fastp/${x}_
    samtools index output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam
done


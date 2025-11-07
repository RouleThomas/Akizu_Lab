#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00

nthreads=12

x=("WT_Rep1" "WT_Rep2" "WT_Rep3" "KO_Rep1" "KO_Rep2" "KO_Rep3")

module load STAR/2.7.3a-GCC-9.3.0
module load SAMtools


# --> UPDATE EACH PATH ACCORDINGLY!!!!!!!!!!!!!!!!!!! 



for x in "${x[@]}"; do
	STAR --genomeDir /Master/meta/STAR_hg38/ \
		--runThreadN 12 \
		--readFilesCommand zcat \
		--readFilesIn output/fastp/${x}_1.fq.gz output/fastp/${x}_2.fq.gz \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix output/STAR/${x}_
    samtools index output/STAR/${x}_Aligned.sortedByCoord.out.bam
done


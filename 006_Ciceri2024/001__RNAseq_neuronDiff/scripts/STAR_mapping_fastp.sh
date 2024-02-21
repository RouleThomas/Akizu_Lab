#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6   # Number of CPU cores per task


nthreads=6

x=(
    "100dN_WT_R1"
    "75dN_WT_R1"
    "50dN_WT_R1"
    "25dN_WT_R1"
    "NPC_WT_R1"
    "ESC_WT_R1"
    "100dN_WT_R2"
    "75dN_WT_R2"
    "50dN_WT_R2"
    "25dN_WT_R2"
    "NPC_WT_R2"
    "ESC_WT_R2"
    "100dN_WT_R3"
    "75dN_WT_R3"
    "50dN_WT_R3"
    "25dN_WT_R3"
    "NPC_WT_R3"
    "ESC_WT_R3"
)

module load STAR/2.7.3a-GCC-9.3.0
module load SAMtools/1.16.1-GCC-11.3.0

for x in "${x[@]}"; do
	STAR --genomeDir ../../Master/meta/STAR_hg38/ \
		--runThreadN 6 \
		--readFilesCommand zcat \
		--readFilesIn output/fastp/${x}_1.fq.gz output/fastp/${x}_2.fq.gz \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix output/STAR/fastp/${x}_
    samtools index output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam
done


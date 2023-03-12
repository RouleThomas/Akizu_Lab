#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G

nthreads=12


module load STAR/2.7.3a-GCC-9.3.0

# example for 1 file:
STAR --genomeDir ../../Master/meta/STAR_hg19/ \
	--runThreadN 12 \
	--readFilesCommand zcat \
	--readFilesIn input/NPC_WT_R1_1.fq.gz input/NPC_WT_R1_2.fq.gz \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix output/STAR/NPC_WT_
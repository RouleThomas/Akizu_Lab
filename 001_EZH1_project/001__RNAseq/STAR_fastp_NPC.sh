#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G
#SBATCH --time=72:00:00

nthreads=12

x=("NPC_WT_R1" "NPC_WT_R2" "NPC_WT_R3"
   "NPC_KO_R1" "NPC_KO_R2" "NPC_KO_R3"
   "NPC_HET_R1" "NPC_HET_R2" "NPC_HET_R3")

module load STAR/2.7.3a-GCC-9.3.0

# example for 1 file:
STAR --genomeDir ../../Master/meta/STAR_hg19/ \
	--runThreadN 12 \
	--readFilesCommand zcat \
	--readFilesIn output/fastp/${x}_1.fq.gz output/fastp/${x}_2.fq.gz \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix output/STAR/fastp/${x}_
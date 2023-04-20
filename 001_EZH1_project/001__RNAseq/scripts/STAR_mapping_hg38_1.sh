#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=250G
#SBATCH --time=100:00:00

nthreads=12

x=("ESC_WT_R1" "ESC_WT_R2" "ESC_WT_R3"
   "ESC_KO_R1" "ESC_KO_R2" "ESC_KO_R3"
   "ESC_HET_R1" "ESC_HET_R2" "ESC_HET_R3"
   "NPC_WT_R1" "NPC_WT_R2" "NPC_WT_R3"
   "NPC_KO_R1" "NPC_KO_R2" "NPC_KO_R3"
   "NPC_HET_R1" "NPC_HET_R2" "NPC_HET_R3"
   "2dN_WT_R1" "2dN_WT_R2" "2dN_WT_R3"
   "2dN_KO_R1" "2dN_KO_R2" "2dN_KO_R3"
   "2dN_HET_R1" "2dN_HET_R2" "2dN_HET_R3")

module load STAR/2.7.3a-GCC-9.3.0
module load sam-bcf-tools/1.6

for x in "${x[@]}"; do
	STAR --genomeDir ../../Master/meta/STAR_hg38/ \
		--runThreadN 12 \
		--readFilesCommand zcat \
		--readFilesIn output/fastp/${x}_1.fq.gz output/fastp/${x}_2.fq.gz \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix output/STAR_hg38/${x}_
    samtools index output/STAR_hg38/${x}_Aligned.sortedByCoord.out.bam
done


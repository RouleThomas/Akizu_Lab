#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=100:00:00

x=("ESC_WT_R1" "ESC_WT_R2" "ESC_WT_R3"
   "ESC_KO_R1" "ESC_KO_R2" "ESC_KO_R3"
   "ESC_HET_R1" "ESC_HET_R2" "ESC_HET_R3"
   "NPC_WT_R1" "NPC_WT_R2" "NPC_WT_R3"
   "NPC_KO_R1" "NPC_KO_R2" "NPC_KO_R3"
   "NPC_HET_R1" "NPC_HET_R2" "NPC_HET_R3"
   "2dN_WT_R1" "2dN_WT_R2" "2dN_WT_R3"
   "2dN_KO_R1" "2dN_KO_R2" "2dN_KO_R3"
   "2dN_HET_R1" "2dN_HET_R2" "2dN_HET_R3")


for x in "${x[@]}"; do
    bamCoverage --bam output/STAR_hg38/${x}_Aligned.sortedByCoord.out.bam \
        --outFileName output/bigwig_hg38/${x}_Aligned.sortedByCoord.out.bw \
        --outFileFormat bigwig \
        --normalizeUsing BPM \
        --binSize 50  
done



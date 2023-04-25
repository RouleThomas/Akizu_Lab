#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=100:00:00

x=("4wN_WT_R1" "4wN_WT_R2" "4wN_KO_R1"
   "4wN_KO_R2" "4wN_HET_R1" "4wN_HET_R2"
   "4wN_HET_R3" "4wN_HET_R4" "4wN_iPSCWT_R1"
   "4wN_iPSCWT_R2" "4wN_iPSCpatient_R1" "4wN_iPSCpatient_R2"
   "8wN_WT_R1" "8wN_WT_R2" "8wN_WT_R3" "8wN_WT_R4" "8wN_KO_R1"
   "8wN_KO_R2" "8wN_KO_R3" "8wN_KO_R4" "8wN_HET_R1" "8wN_HET_R2"
   "8wN_HET_R3" "8wN_HET_R4" "8wN_iPSCpatient_R1" "8wN_iPSCpatient_R2")


for x in "${x[@]}"; do
    bamCoverage --bam output/STAR_hg38/${x}_Aligned.sortedByCoord.out.bam \
        --outFileName output/bigwig_hg38/${x}_Aligned.sortedByCoord.out.bw \
        --outFileFormat bigwig \
        --normalizeUsing BPM \
        --binSize 50  
done



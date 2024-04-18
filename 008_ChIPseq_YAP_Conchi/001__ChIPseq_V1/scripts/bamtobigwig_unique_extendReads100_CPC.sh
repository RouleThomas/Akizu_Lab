#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=(
"CPC_untreated_YAP1_R1"
"CPC_RA_YAP1_R1"
"CPC_untreated_YAP1_R2"
"CPC_RA_YAP1_R2"
"CPC_untreated_TEAD4_R1"
"CPC_RA_TEAD4_R1"
"CPC_untreated_TEAD4_R2"
"CPC_RA_TEAD4_R2"
"CPC_untreated_NR2F2_R1"
"CPC_RA_NR2F2_R1"
"CPC_untreated_NR2F2_R2"
"CPC_RA_NR2F2_R2"
"CPC_RA_NR2F2_R3"
"CPC_untreated_NR2F2_R3"
"CPC_untreated_input_R3"
"CPC_RA_input_R3"
"CPC_untreated_YAP1_R3"
"CPC_RA_YAP1_R3"
"CPC_untreated_TEAD4_R3"
"CPC_RA_TEAD4_R3"
)

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.unique.dupmark.sorted.bam \
        --outFileName output/bigwig_noextendReads/${x}.unique.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads 100 \
        --scaleFactor 0.5
done


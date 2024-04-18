#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


bamCoverage --bam output/bowtie2/xxx.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/xxx.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads xxx \
    --scaleFactor 0.5


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




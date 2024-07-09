#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=(
"hESC_WT_input_R4"
"hESC_WT_input_R5"
"hESC_YAPKO_input_R4"
"hESC_WT_QSER1_R4"
"hESC_WT_QSER1_R5"
"hESC_YAPKO_QSER1_R4"
"hESC_YAPKO_QSER1_R5"
"hESC_WT_TEAD4_R4"
"hESC_WT_TEAD4_R5"
"hESC_WT_YAP1_R4"
"hESC_WT_YAP1_R5"
)

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.unique.dupmark.sorted.bam \
        --outFileName output/bigwig_extendReads100/${x}.unique.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads 100 \
        --scaleFactor 0.5
done


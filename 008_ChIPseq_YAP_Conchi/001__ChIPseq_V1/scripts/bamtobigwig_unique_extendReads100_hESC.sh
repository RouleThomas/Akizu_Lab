#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=(
"hESC_WT_DVL2_R1"
"hESC_YAPKO_DVL2_R1"
"hESC_WT_EZH2_R1"
"hESC_YAPKO_EZH2_R1"
"hESC_WT_EZH2_R2"
"hESC_YAPKO_EZH2_R2"
"hESC_WT_QSER1_R1"
"hESC_YAPKO_QSER1_R1"
"hESC_WT_QSER1_R2"
"hESC_YAPKO_QSER1_R2"
"hESC_WT_input_R1"
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


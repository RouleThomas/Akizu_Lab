#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7


input_list=(
    "PSC_WT_EZH1_006R"
    "PSC_WT_EZH2_006R"
    "PSC_WT_EZH2_010R"
    "PSC_WT_EZH2_014R1"
    "PSC_WT_H3K27me3_006R"
    "PSC_WT_H3K27me3_010R"
    "PSC_WT_H3K27me3_013R1"
    "PSC_WT_IGG_006R"
    "PSC_WT_IGG_010R"
    "PSC_WT_IGG_013R1"
    "PSC_WT_IGG_014R1"
    "PSC_WT_SUZ12_006R"
    "PSC_WT_SUZ12_013R1"
    "PSC_WT_SUZ12_014R1"
)

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.unique.dupmark.sorted.bam \
        --outFileName output/bigwig/${x}.unique.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done


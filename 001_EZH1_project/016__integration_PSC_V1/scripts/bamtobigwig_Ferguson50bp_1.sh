#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7


input_list=(
    "PSC_KOEF1aEZH1_EZH1_005R"
    "PSC_KOEF1aEZH1_EZH1_006R"
    "PSC_KOEF1aEZH1_EZH1_013R1"
    "PSC_KOEF1aEZH1_EZH2_006R"
    "PSC_KOEF1aEZH1_EZH2_013R1"
    "PSC_KOEF1aEZH1_EZH2_014R1"
    "PSC_KOEF1aEZH1_H3K27me3_005R"
    "PSC_KOEF1aEZH1_H3K27me3_006R"
    "PSC_KOEF1aEZH1_H3K27me3_013R1"
    "PSC_KOEF1aEZH1_IGG_005R"
    "PSC_KOEF1aEZH1_IGG_006R"
    "PSC_KOEF1aEZH1_IGG_013R1"
    "PSC_KOEF1aEZH1_IGG_014R1"
    "PSC_KOEF1aEZH1_SUZ12_005R"
    "PSC_KOEF1aEZH1_SUZ12_006R"
    "PSC_KOEF1aEZH1_SUZ12_013R1"
)

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.unique.dupmark.sorted.bam \
        --outFileName output/bigwig_Ferguson50bp/${x}.50bp.unique.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 50 \
        --numberOfProcessors 7 \
        --extendReads
done


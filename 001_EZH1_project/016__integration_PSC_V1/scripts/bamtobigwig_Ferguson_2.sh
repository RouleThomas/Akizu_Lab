#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7


input_list=(
    "PSC_KO_EZH1_006R"
    "PSC_KO_EZH1_013R1"
    "PSC_KO_EZH1_014R2"
    "PSC_KO_EZH2_013R1"
    "PSC_KO_EZH2_014R1"
    "PSC_KO_EZH2_014R2"
    "PSC_KO_H3K27me3_006R"
    "PSC_KO_H3K27me3_013R1"
    "PSC_KO_H3K27me3_014R2"
    "PSC_KO_IGG_006R"
    "PSC_KO_IGG_013R1"
    "PSC_KO_IGG_014R1"
    "PSC_KO_IGG_014R2"
    "PSC_KO_SUZ12_013R1"
    "PSC_KO_SUZ12_014R1"
    "PSC_KO_SUZ12_014R2"
)

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.filter.bam \
        --outFileName output/bigwig/${x}.filter.bw \
        --outFileFormat bigwig \
        --binSize 50 \
        --numberOfProcessors 7 \
        --extendReads
done


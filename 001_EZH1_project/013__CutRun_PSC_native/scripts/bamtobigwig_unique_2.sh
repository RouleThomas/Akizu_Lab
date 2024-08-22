#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=(
"PSC_WT_IGG_R1"
"PSC_WT_SUZ12_R1"
"PSC_KO_EZH1_R2"
"PSC_KO_EZH2_R2"
"PSC_KO_H3K27me3_R2"
"PSC_KO_IGG_R2"
"PSC_KO_SUZ12_R2"
"PSC_KOEF1aEZH1_EZH1_R2"
"PSC_KOEF1aEZH1_EZH2_R2"
"PSC_KOEF1aEZH1_H3K27me3_R2"
"PSC_KOEF1aEZH1_IGG_R2"
"PSC_KOEF1aEZH1_SUZ12_R2"
"PSC_WTEF1aEZH1_EZH1_R2"
"PSC_WTEF1aEZH1_EZH2_R2"
"PSC_WTEF1aEZH1_H3K27me3_R2"
"PSC_WTEF1aEZH1_IGG_R2"
"PSC_WTEF1aEZH1_SUZ12_R2"
"PSC_WT_EZH1_R2"
"PSC_WT_EZH2_R2"
"PSC_WT_H3K27me3_R2"
"PSC_WT_IGG_R2"
"PSC_WT_SUZ12_R2"
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

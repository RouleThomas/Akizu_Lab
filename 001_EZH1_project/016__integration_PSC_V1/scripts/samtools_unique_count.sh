#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00

nthreads=5

module load SAMtools/1.16.1-GCC-11.3.0
module load picard/2.26.10-Java-15

x=(
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


for x in "${x[@]}"; do
    # count the uniq map reads
    samtools view -S -F 4 -c output/bowtie2/${x}.unique.dupmark.sorted.bam > output/bowtie2/${x}.unique.dupmark.sorted-count.txt
done



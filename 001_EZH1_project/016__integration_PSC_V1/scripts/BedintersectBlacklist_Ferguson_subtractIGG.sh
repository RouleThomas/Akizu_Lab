#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00


input_list=(

    "PSC_KOEF1aEZH1_EZH2_006R"
    "PSC_KOEF1aEZH1_EZH2_013R1"
    "PSC_KOEF1aEZH1_EZH2_014R1"
    "PSC_KOEF1aEZH1_H3K27me3_005R"
    "PSC_KOEF1aEZH1_H3K27me3_006R"
    "PSC_KOEF1aEZH1_H3K27me3_013R1"
    "PSC_KOEF1aEZH1_SUZ12_005R"
    "PSC_KOEF1aEZH1_SUZ12_006R"
    "PSC_KOEF1aEZH1_SUZ12_013R1"


    "PSC_KO_EZH2_013R1"
    "PSC_KO_EZH2_014R1"
    "PSC_KO_EZH2_014R2"
    "PSC_KO_H3K27me3_006R"
    "PSC_KO_H3K27me3_013R1"
    "PSC_KO_H3K27me3_014R2"
    "PSC_KO_SUZ12_013R1"
    "PSC_KO_SUZ12_014R1"
    "PSC_KO_SUZ12_014R2"

    "PSC_WT_EZH2_006R"
    "PSC_WT_EZH2_010R"
    "PSC_WT_EZH2_014R1"
    "PSC_WT_H3K27me3_006R"
    "PSC_WT_H3K27me3_010R"
    "PSC_WT_H3K27me3_013R1"
    "PSC_WT_SUZ12_006R"
    "PSC_WT_SUZ12_013R1"
    "PSC_WT_SUZ12_014R1"
)


for x in "${input_list[@]}"; do
    bedtools intersect -v \
        -a output/bigwig/${x}_subtract.filter.bedGraph \
        -b ../../Master/meta/hg38-blacklist.v2.bed > output/bigwig/${x}_subtract.filter.blacklist.bedGraph
done




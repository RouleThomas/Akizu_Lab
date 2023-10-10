#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=25G


macs2_out="output/macs2"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"

# average basepair q-value threshold (log5)
q=5

# make output directory
mkdir -p ${macs2_out}/broad_noMask_qval${q}
mkdir -p ${macs2_out}/broad_blacklist_qval${q}

  

x=("NPC_WT_EZH1cs"
    "NPC_WT_EZH1pt"
    "NPC_WT_EZH2"
    "NPC_WT_H3K27me1"
    "NPC_WT_H3K27me3" 
    "NPC_WT_SUZ12"
    "NPC_KO_EZH1cs" 
    "NPC_KO_EZH1pt"
    "NPC_KO_EZH2"
    "NPC_KO_H3K27me1" 
    "NPC_KO_H3K27me3" 
    "NPC_KO_SUZ12"

    "PSC_KOEF1aEZH1_EZH1cs" 
    "PSC_KOEF1aEZH1_EZH1pt" 
    "PSC_KOEF1aEZH1_H3K27me3" 
    "PSC_KOEF1aEZH1_HA" 
    "PSC_KOEF1aEZH1_SUZ12" 
    "PSC_KOsynEZH1_EZH1cs" 
    "PSC_KOsynEZH1_EZH1pt" 
    "PSC_KOsynEZH1_H3K27me3" 
    "PSC_KOsynEZH1_HA" 
    "PSC_KOsynEZH1_SUZ12")



for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.broadPeak > ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak \
        -b ${blacklist} > ${macs2_out}/broad_blacklist_qval${q}/${x}_peaks.broadPeak   
done








q=5

x=("NPC_WT_H3K4me3"
"NPC_KO_H3K4me3")


for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.broadPeak > ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak \
        -b ${blacklist} > ${macs2_out}/broad_blacklist_qval${q}/${x}_peaks.broadPeak   
done
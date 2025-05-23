#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"

# average basepair q-value threshold (log5)
q=5  #1.30103

# make output directory
mkdir -p ${macs2_out}/broad_noMask_qval${q}
mkdir -p ${macs2_out}/broad_blacklist_qval${q}

  

x=(
    "PSC_KOEF1aEZH1_EZH1_pool"
    "PSC_KOEF1aEZH1_EZH2_pool"
    "PSC_KOEF1aEZH1_H3K27me3_pool"
    "PSC_KOEF1aEZH1_SUZ12_pool"

    "PSC_KO_EZH1_pool"
    "PSC_KO_EZH2_pool"
    "PSC_KO_H3K27me3_pool"
    "PSC_KO_SUZ12_pool"

    "PSC_WT_EZH2_pool"
    "PSC_WT_H3K27me3_pool"
    "PSC_WT_SUZ12_pool"
)



for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.broadPeak > ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak \
        -b ${blacklist} > ${macs2_out}/broad_blacklist_qval${q}/${x}_peaks.broadPeak   
done

#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=25G


macs2_out="output/macs2"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"

# average basepair q-value threshold (log5)
q=1.30103

# make output directory
mkdir -p ${macs2_out}/broad_noMask_qval${q}
mkdir -p ${macs2_out}/broad_blacklist_qval${q}


x=("8wN_HET_H3K27me3_R1"
    "8wN_HET_H3K27me3_R2" 
    "8wN_HET_H3K27me3_R3" 
    "8wN_HET_H3K27me3_R4" 
    "8wN_KO_H3K27me3_R1" 
    "8wN_KO_H3K27me3_R2" 
    "8wN_KO_H3K27me3_R3" 
    "8wN_KO_H3K27me3_R4" 
    "8wN_WT_H3K27me3_R1" 
    "8wN_WT_H3K27me3_R2"
    "8wN_WT_H3K27me3_R3" 
    "8wN_WT_H3K27me3_R4" 
    "8wN_iPSCpatient_H3K27me3_R1" 
    "8wN_iPSCpatient_H3K27me3_R2")



for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.broadPeak > ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak \
        -b ${blacklist} > ${macs2_out}/broad_blacklist_qval${q}/${x}_peaks.broadPeak   
done

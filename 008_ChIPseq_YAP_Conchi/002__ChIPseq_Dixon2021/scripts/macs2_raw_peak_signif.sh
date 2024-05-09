#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=25G



## broad pool


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"

# average basepair q-value threshold (log5)
q=1.30103 # 1.30103

# make output directory
mkdir -p ${macs2_out}/broad_noMask_qval${q}
mkdir -p ${macs2_out}/broad_blacklist_qval${q}



x=( 
    "hESC_WT_QSER1FLAG_pool"
    "hESC_WT_DNMT3A_pool"
    "hESC_WT_DNMT3B_pool"
    "hESC_WT_H3K27me3_pool"
    "hESC_WT_H3K4me3_pool"
 )




for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.broadPeak > ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak \
        -b ${blacklist} > ${macs2_out}/broad_blacklist_qval${q}/${x}_peaks.broadPeak   
done





## broad replicate 


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"

# average basepair q-value threshold (log5)
q=4 # 1.30103

# make output directory
mkdir -p ${macs2_out}/broad_noMask_qval${q}
mkdir -p ${macs2_out}/broad_blacklist_qval${q}




x=( 
    "hESC_WT_QSER1FLAG_R1"
    "hESC_WT_QSER1FLAG_R2"
    "hESC_WT_DNMT3A_R1"
    "hESC_WT_DNMT3A_R2"
    "hESC_WT_DNMT3B_R1"
    "hESC_WT_DNMT3B_R2"
    "hESC_WT_H3K27me3_R1"
    "hESC_WT_H3K27me3_R2"
    "hESC_WT_H3K4me3_R1"
    "hESC_WT_H3K4me3_R2"
    "hESC_WT_EZH2_R1"
 )


for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.broadPeak > ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak \
        -b ${blacklist} > ${macs2_out}/broad_blacklist_qval${q}/${x}_peaks.broadPeak   
done





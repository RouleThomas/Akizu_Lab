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


x=("2dN_HET_H3K27me3_pool"
"2dN_KO_H3K27me3_pool"
"2dN_WT_H3K27me3_pool"
"ESC_HET_H3K27me3_pool"
"ESC_KO_H3K27me3_pool"
"ESC_WT_H3K27me3_pool"
"NPC_HET_H3K27me3_pool"
"NPC_KO_H3K27me3_pool"
"NPC_WT_H3K27me3_pool")



for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.broadPeak > ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak \
        -b ${blacklist} > ${macs2_out}/broad_blacklist_qval${q}/${x}_peaks.broadPeak   
done

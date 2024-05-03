#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=25G




## broad replicate 


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"

# average basepair q-value threshold (log5)
q=5 # 1.30103

# make output directory
mkdir -p ${macs2_out}/broad_noMask_qval${q}
mkdir -p ${macs2_out}/broad_blacklist_qval${q}


x=( 
    "hESC_WT_TEAD4_R1"
    "hESC_WT_YAP1_R1"
 )


for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.broadPeak > ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak \
        -b ${blacklist} > ${macs2_out}/broad_blacklist_qval${q}/${x}_peaks.broadPeak   
done





## narrow replicate 

macs2_out="output/macs2/narrow"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"

# average basepair q-value threshold (log5)
q=5 # 1.30103

# make output directory
mkdir -p ${macs2_out}/narrow_noMask_qval${q}
mkdir -p ${macs2_out}/narrow_blacklist_qval${q}

x=( 
    "hESC_WT_TEAD4_R1"
    "hESC_WT_YAP1_R1"
 )


for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.narrowPeak > ${macs2_out}/narrow_noMask_qval${q}/${x}_peaks.narrowPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/narrow_noMask_qval${q}/${x}_peaks.narrowPeak \
        -b ${blacklist} > ${macs2_out}/narrow_blacklist_qval${q}/${x}_peaks.narrowPeak   
done






















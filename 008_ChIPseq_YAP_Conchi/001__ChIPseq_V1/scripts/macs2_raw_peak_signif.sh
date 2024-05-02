#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=25G



## broad pool


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"

# average basepair q-value threshold (log5)
q=4 # 1.30103

# make output directory
mkdir -p ${macs2_out}/broad_noMask_qval${q}
mkdir -p ${macs2_out}/broad_blacklist_qval${q}



x=( 
    "hESC_WT_EZH2_pool"
    "hESC_YAPKO_EZH2_pool"
    "hESC_WT_QSER1_pool"
    "hESC_YAPKO_QSER1_pool"

    "CPC_untreated_NR2F2_pool"
    "CPC_RA_NR2F2_pool"
 )




for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.broadPeak > ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak \
        -b ${blacklist} > ${macs2_out}/broad_blacklist_qval${q}/${x}_peaks.broadPeak   
done



## narrow pool


macs2_out="output/macs2/narrow"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"

# average basepair q-value threshold (log5)
q=1.30103 # 1.30103

# make output directory
mkdir -p ${macs2_out}/narrow_noMask_qval${q}
mkdir -p ${macs2_out}/narrow_blacklist_qval${q}



x=( 
    "hESC_WT_EZH2_pool"
    "hESC_YAPKO_EZH2_pool"
    "hESC_WT_QSER1_pool"
    "hESC_YAPKO_QSER1_pool"

    "CPC_untreated_NR2F2_pool"
    "CPC_RA_NR2F2_pool"
 )



for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.narrowPeak > ${macs2_out}/narrow_noMask_qval${q}/${x}_peaks.narrowPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/narrow_noMask_qval${q}/${x}_peaks.narrowPeak \
        -b ${blacklist} > ${macs2_out}/narrow_blacklist_qval${q}/${x}_peaks.narrowPeak   
done





## broad replicate 


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"

# average basepair q-value threshold (log5)
q=2.30103 # 1.30103

# make output directory
mkdir -p ${macs2_out}/broad_noMask_qval${q}
mkdir -p ${macs2_out}/broad_blacklist_qval${q}




x=( 
    "hESC_WT_DVL2_R1"
    "hESC_YAPKO_DVL2_R1"
    "hESC_WT_EZH2_R1"
    "hESC_YAPKO_EZH2_R1"
    "hESC_WT_EZH2_R2"
    "hESC_YAPKO_EZH2_R2"
    "hESC_WT_QSER1_R1"
    "hESC_YAPKO_QSER1_R1"
    "hESC_WT_QSER1_R2"
    "hESC_YAPKO_QSER1_R2"

    "CPC_untreated_YAP1_R1"
    "CPC_untreated_YAP1_R2"
    "CPC_untreated_TEAD4_R1"
    "CPC_untreated_TEAD4_R2"
    "CPC_untreated_NR2F2_R1"
    "CPC_untreated_NR2F2_R2"
    "CPC_untreated_NR2F2_R3"
    "CPC_untreated_YAP1_R3"
    "CPC_untreated_TEAD4_R3"
    "CPC_RA_YAP1_R1"
    "CPC_RA_YAP1_R2"
    "CPC_RA_TEAD4_R1"
    "CPC_RA_TEAD4_R2"
    "CPC_RA_NR2F2_R1"
    "CPC_RA_NR2F2_R2"
    "CPC_RA_NR2F2_R3"

    "CPC_RA_TEAD4_R3"
    "CPC_RA_YAP1_R3"
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
q=3 # 1.30103

# make output directory
mkdir -p ${macs2_out}/narrow_noMask_qval${q}
mkdir -p ${macs2_out}/narrow_blacklist_qval${q}




x=( 
    "hESC_WT_DVL2_R1"
    "hESC_YAPKO_DVL2_R1"
    "hESC_WT_EZH2_R1"
    "hESC_YAPKO_EZH2_R1"
    "hESC_WT_EZH2_R2"
    "hESC_YAPKO_EZH2_R2"
    "hESC_WT_QSER1_R1"
    "hESC_YAPKO_QSER1_R1"
    "hESC_WT_QSER1_R2"
    "hESC_YAPKO_QSER1_R2"

    "CPC_untreated_YAP1_R1"
    "CPC_untreated_YAP1_R2"
    "CPC_untreated_TEAD4_R1"
    "CPC_untreated_TEAD4_R2"
    "CPC_untreated_NR2F2_R1"
    "CPC_untreated_NR2F2_R2"
    "CPC_untreated_NR2F2_R3"
    "CPC_untreated_YAP1_R3"
    "CPC_untreated_TEAD4_R3"
    "CPC_RA_YAP1_R1"
    "CPC_RA_YAP1_R2"
    "CPC_RA_TEAD4_R1"
    "CPC_RA_TEAD4_R2"
    "CPC_RA_NR2F2_R1"
    "CPC_RA_NR2F2_R2"
    "CPC_RA_NR2F2_R3"

    "CPC_RA_TEAD4_R3"
    "CPC_RA_YAP1_R3"
 )



for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.narrowPeak > ${macs2_out}/narrow_noMask_qval${q}/${x}_peaks.narrowPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/narrow_noMask_qval${q}/${x}_peaks.narrowPeak \
        -b ${blacklist} > ${macs2_out}/narrow_blacklist_qval${q}/${x}_peaks.narrowPeak   
done




























macs2_out="output/macs2/broad"

macs2_out="output/macs2/narrow"

blacklist="../../Master/meta/hg38-blacklist.v2.bed"


# average basepair q-value threshold (log5)
q=1.30103 # 1.30103

# make output directory
mkdir -p ${macs2_out}/broad_noMask_qval${q}
mkdir -p ${macs2_out}/broad_blacklist_qval${q}

mkdir -p ${macs2_out}/narrow_noMask_qval${q}
mkdir -p ${macs2_out}/narrow_blacklist_qval${q}


x=( 
    "hESC_WT_DVL2_R1"
    "hESC_YAPKO_DVL2_R1"
    "hESC_WT_EZH2_R1"
    "hESC_YAPKO_EZH2_R1"
    "hESC_WT_EZH2_R2"
    "hESC_YAPKO_EZH2_R2"
    "hESC_WT_QSER1_R1"
    "hESC_YAPKO_QSER1_R1"
    "hESC_WT_QSER1_R2"
    "hESC_YAPKO_QSER1_R2"

    "CPC_untreated_YAP1_R1"
    "CPC_untreated_YAP1_R2"
    "CPC_untreated_TEAD4_R1"
    "CPC_untreated_TEAD4_R2"
    "CPC_untreated_NR2F2_R1"
    "CPC_untreated_NR2F2_R2"
    "CPC_untreated_NR2F2_R3"
    "CPC_untreated_YAP1_R3"
    "CPC_untreated_TEAD4_R3"
    "CPC_RA_YAP1_R1"
    "CPC_RA_YAP1_R2"
    "CPC_RA_TEAD4_R1"
    "CPC_RA_TEAD4_R2"
    "CPC_RA_NR2F2_R1"
    "CPC_RA_NR2F2_R2"
    "CPC_RA_NR2F2_R3"


 )

    "CPC_RA_TEAD4_R3"
    "CPC_RA_YAP1_R3"

# pool samples
x=( 
    "hESC_WT_EZH2_pool"
    "hESC_YAPKO_EZH2_pool"


    "hESC_WT_QSER1_pool"
    "hESC_YAPKO_QSER1_pool"


    "CPC_untreated_NR2F2_pool"
    "CPC_RA_NR2F2_pool"

 )



for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.narrowPeak > ${macs2_out}/narrow_noMask_qval${q}/${x}_peaks.narrowPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/narrow_noMask_qval${q}/${x}_peaks.narrowPeak \
        -b ${blacklist} > ${macs2_out}/narrow_blacklist_qval${q}/${x}_peaks.narrowPeak   
done




for x in "${x[@]}"; do
    awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q {print}' ${macs2_out}/${x}_peaks.broadPeak > ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak                                             
    
    bedtools intersect -v -wa \
        -a ${macs2_out}/broad_noMask_qval${q}/${x}_peaks.broadPeak \
        -b ${blacklist} > ${macs2_out}/broad_blacklist_qval${q}/${x}_peaks.broadPeak   
done






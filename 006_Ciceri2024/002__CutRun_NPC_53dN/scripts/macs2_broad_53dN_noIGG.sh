#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad_noIGG"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"




input_list=(
    "53dN_WT_H3K27ac_R1"
    "53dN_WT_H3K27ac_R2"
    "53dN_WT_H3K27me3_R1"
    "53dN_WT_H3K27me3_R2"
    "53dN_WT_H3K4me3_R1"
    "53dN_WT_H3K4me3_R2"
    "53dN_WT_H3K9me3_R1"
    "53dN_WT_H3K9me3_R2"
    "53dN_WT_IGG_R1"
    "53dN_WT_IGG_R2"
)

for x in "${input_list[@]}"; do
    macs2 callpeak -t output/bowtie2/${x}.unique.dupmark.sorted.bam \
            -f BAMPE --keep-dup auto \
            --nomodel -g hs \
            --outdir ${macs2_out} -n ${x} --broad 

done


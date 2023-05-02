#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=150G


macs2_out="output/macs2"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


macs2 callpeak -t output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.bam \
    output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.bam \
    output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.bam \
    output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.bam \
    -c output/bowtie2/8wN_KO_IGG_R1.dupmark.sorted.bam \
    output/bowtie2/8wN_KO_IGG_R2.dupmark.sorted.bam \
    output/bowtie2/8wN_KO_IGG_R3.dupmark.sorted.bam \
    output/bowtie2/8wN_KO_IGG_R4.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n 8wN_KO_H3K27me3_pool --broad 


macs2 callpeak -t output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.bam \
    output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.bam \
    output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.bam \
    output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.bam \
    -c output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.bam \
    output/bowtie2/8wN_HET_IGG_R2.dupmark.sorted.bam \
    output/bowtie2/8wN_HET_IGG_R3.dupmark.sorted.bam \
    output/bowtie2/8wN_HET_IGG_R4.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n 8wN_HET_H3K27me3_pool --broad 

macs2 callpeak -t output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.bam \
    output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.bam \
    output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.bam \
    output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.bam \
    -c output/bowtie2/8wN_WT_IGG_R1.dupmark.sorted.bam \
    output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.bam \
    output/bowtie2/8wN_WT_IGG_R3.dupmark.sorted.bam \
    output/bowtie2/8wN_WT_IGG_R4.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n 8wN_WT_H3K27me3_pool --broad 


macs2 callpeak -t output/bowtie2/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.bam \
    -c output/bowtie2/8wN_iPSCpatient_IGG_R1.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n 8wN_iPSCpatient_H3K27me3_pool --broad 



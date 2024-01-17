#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"


macs2 callpeak -t output/bowtie2/50dN_WTQ731E_H3K27me3_R1.unique.dupmark.sorted.bam \
    output/bowtie2/50dN_WTQ731E_H3K27me3_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/50dN_WTQ731E_IGG_R1.unique.dupmark.sorted.bam \
    output/bowtie2/50dN_WTQ731E_IGG_R2.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n 50dN_WTQ731E_H3K27me3_pool --broad 


macs2 callpeak -t output/bowtie2/50dN_KO_H3K27me3_R1.unique.dupmark.sorted.bam \
    output/bowtie2/50dN_KO_H3K27me3_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/50dN_KO_IGG_R1.unique.dupmark.sorted.bam \
    output/bowtie2/50dN_KO_IGG_R2.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n 50dN_KO_H3K27me3_pool --broad 

macs2 callpeak -t output/bowtie2/50dN_KOEF1aEZH1_H3K27me3_R1.unique.dupmark.sorted.bam \
    output/bowtie2/50dN_KOEF1aEZH1_H3K27me3_R2.unique.dupmark.sorted.bam \
    -c output/bowtie2/50dN_KOEF1aEZH1_IGG_R1.unique.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n 50dN_KOEF1aEZH1_H3K27me3_pool --broad 


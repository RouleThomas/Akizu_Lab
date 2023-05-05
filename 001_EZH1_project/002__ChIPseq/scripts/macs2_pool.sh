#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=150G


macs2_out="output/macs2"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"



macs2 callpeak -t output/bowtie2_endtoend/2dN_HET_H3K27me3_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/2dN_HET_H3K27me3_R2.dupmark.sorted.bam \
    -c output/bowtie2_endtoend/2dN_HET_input_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/2dN_HET_input_R2.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n 2dN_HET_H3K27me3_pool --broad 


macs2 callpeak -t output/bowtie2_endtoend/2dN_KO_H3K27me3_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/2dN_KO_H3K27me3_R2.dupmark.sorted.bam \
    -c output/bowtie2_endtoend/2dN_KO_input_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/2dN_KO_input_R2.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n 2dN_KO_H3K27me3_pool --broad 


macs2 callpeak -t output/bowtie2_endtoend/2dN_WT_H3K27me3_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/2dN_WT_H3K27me3_R2.dupmark.sorted.bam \
    -c output/bowtie2_endtoend/2dN_WT_input_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/2dN_WT_input_R2.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n 2dN_WT_H3K27me3_pool --broad 


macs2 callpeak -t output/bowtie2_endtoend/ESC_HET_H3K27me3_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/ESC_HET_H3K27me3_R2.dupmark.sorted.bam \
    -c output/bowtie2_endtoend/ESC_HET_input_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/ESC_HET_input_R2.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n ESC_HET_H3K27me3_pool --broad 

macs2 callpeak -t output/bowtie2_endtoend/ESC_KO_H3K27me3_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/ESC_KO_H3K27me3_R2.dupmark.sorted.bam \
    -c output/bowtie2_endtoend/ESC_KO_input_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/ESC_KO_input_R2.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n ESC_KO_H3K27me3_pool --broad 

macs2 callpeak -t output/bowtie2_endtoend/ESC_WT_H3K27me3_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/ESC_WT_H3K27me3_R2.dupmark.sorted.bam \
    -c output/bowtie2_endtoend/ESC_WT_input_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/ESC_WT_input_R2.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n ESC_WT_H3K27me3_pool --broad 


macs2 callpeak -t output/bowtie2_endtoend/NPC_HET_H3K27me3_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/NPC_HET_H3K27me3_R2.dupmark.sorted.bam \
    -c output/bowtie2_endtoend/NPC_HET_input_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/NPC_HET_input_R2.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n NPC_HET_H3K27me3_pool --broad 


macs2 callpeak -t output/bowtie2_endtoend/NPC_KO_H3K27me3_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/NPC_KO_H3K27me3_R2.dupmark.sorted.bam \
    -c output/bowtie2_endtoend/NPC_KO_input_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/NPC_KO_input_R2.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n NPC_KO_H3K27me3_pool --broad 


macs2 callpeak -t output/bowtie2_endtoend/NPC_WT_H3K27me3_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/NPC_WT_H3K27me3_R2.dupmark.sorted.bam \
    -c output/bowtie2_endtoend/NPC_WT_input_R1.dupmark.sorted.bam \
    output/bowtie2_endtoend/NPC_WT_input_R2.dupmark.sorted.bam \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir ${macs2_out} -n NPC_WT_H3K27me3_pool --broad 





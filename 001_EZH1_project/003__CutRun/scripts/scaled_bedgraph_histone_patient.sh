#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=10G



awk -v scaling_factor=2.26817010309278 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=1.47564432989691 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.fragments.histonescaled.bedgraph


awk -v scaling_factor=1.60606060606061 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_iPSCpatient_IGG_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_iPSCpatient_IGG_R1.dupmark.sorted.fragments.histonescaled.bedgraph

awk -v scaling_factor=1.18050065876153 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_iPSCpatient_IGG_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_iPSCpatient_IGG_R2.dupmark.sorted.fragments.histonescaled.bedgraph

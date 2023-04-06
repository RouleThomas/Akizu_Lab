#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=10G



awk -v scaling_factor=1.58479381443299 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=1.88492268041237 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=1.05090206185567 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.fragments.histonescaled.bedgraph




awk -v scaling_factor=2.85734536082474 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=1.38471673254282 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=1.13570487483531 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_IGG_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_IGG_R2.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=1.05270092226614 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_IGG_R3.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_IGG_R3.dupmark.sorted.fragments.histonescaled.bedgraph


awk -v scaling_factor=2.55994729907773 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_IGG_R4.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_IGG_R4.dupmark.sorted.fragments.histonescaled.bedgraph


awk -v scaling_factor=1.30657216494845 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=1 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=3.80953608247423 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.fragments.histonescaled.bedgraph




awk -v scaling_factor=1.12899484536082 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.fragments.histonescaled.bedgraph





awk -v scaling_factor=1.3768115942029 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_IGG_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_IGG_R1.dupmark.sorted.fragments.histonescaled.bedgraph





awk -v scaling_factor=1 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_IGG_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_IGG_R2.dupmark.sorted.fragments.histonescaled.bedgraph


awk -v scaling_factor=1.42555994729908 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_IGG_R3.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_IGG_R3.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=1.02766798418972 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_IGG_R4.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_IGG_R4.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=1.69085051546392 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.fragments.histonescaled.bedgraph


awk -v scaling_factor=2.92551546391753 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=2.05154639175258 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=3.09368556701031 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.fragments.histonescaled.bedgraph


awk -v scaling_factor=2.1870882740448 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_IGG_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_IGG_R1.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=2.80237154150198 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.fragments.histonescaled.bedgraph


awk -v scaling_factor=1.20026350461133 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_IGG_R3.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_IGG_R3.dupmark.sorted.fragments.histonescaled.bedgraph



awk -v scaling_factor=5.21607378129117 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_IGG_R4.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_IGG_R4.dupmark.sorted.fragments.histonescaled.bedgraph

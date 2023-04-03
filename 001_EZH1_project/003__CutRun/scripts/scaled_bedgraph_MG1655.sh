#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=10G



awk -v scaling_factor=2.07929383602633 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.fragments.MG1655scaled.bedgraph

awk -v scaling_factor=1.92654099341712 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.fragments.MG1655scaled.bedgraph

awk -v scaling_factor=1.18237582286056 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.fragments.MG1655scaled.bedgraph


awk -v scaling_factor=2.090065828845 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.fragments.MG1655scaled.bedgraph





awk -v scaling_factor=2.0094484954298 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.fragments.MG1655scaled.bedgraph




awk -v scaling_factor=2.39088014788949 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_IGG_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_IGG_R2.dupmark.sorted.fragments.MG1655scaled.bedgraph


awk -v scaling_factor=1.15040566909726 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_IGG_R3.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_IGG_R3.dupmark.sorted.fragments.MG1655scaled.bedgraph



awk -v scaling_factor=1 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_HET_IGG_R4.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_HET_IGG_R4.dupmark.sorted.fragments.MG1655scaled.bedgraph



awk -v scaling_factor=3.22396768402154 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.fragments.MG1655scaled.bedgraph



awk -v scaling_factor=1.44389587073609 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.fragments.MG1655scaled.bedgraph



awk -v scaling_factor=7.20436864153202 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.fragments.MG1655scaled.bedgraph



awk -v scaling_factor=1 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.fragments.MG1655scaled.bedgraph



awk -v scaling_factor=3.77749820273185 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_IGG_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_IGG_R1.dupmark.sorted.fragments.MG1655scaled.bedgraph




awk -v scaling_factor=1.99388928828181 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_IGG_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_IGG_R2.dupmark.sorted.fragments.MG1655scaled.bedgraph




awk -v scaling_factor=3.34148094895758 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_IGG_R3.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_IGG_R3.dupmark.sorted.fragments.MG1655scaled.bedgraph



awk -v scaling_factor=1.22712334394577 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_KO_IGG_R4.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_KO_IGG_R4.dupmark.sorted.fragments.MG1655scaled.bedgraph



awk -v scaling_factor=2.99640933572711 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.fragments.MG1655scaled.bedgraph


awk -v scaling_factor=5.95347097546379 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.fragments.MG1655scaled.bedgraph



awk -v scaling_factor=3.02154398563734 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.fragments.MG1655scaled.bedgraph





awk -v scaling_factor=6.26256732495512 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.fragments.MG1655scaled.bedgraph





awk -v scaling_factor=3.26091198521105 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_IGG_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_IGG_R1.dupmark.sorted.fragments.MG1655scaled.bedgraph






awk -v scaling_factor=7.21351545650611 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.fragments.MG1655scaled.bedgraph






awk -v scaling_factor=1.8490294751977 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_IGG_R3.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_IGG_R3.dupmark.sorted.fragments.MG1655scaled.bedgraph






awk -v scaling_factor=9.98243812262504 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_WT_IGG_R4.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_WT_IGG_R4.dupmark.sorted.fragments.MG1655scaled.bedgraph




awk -v scaling_factor=2.0050867743866 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.fragments.MG1655scaled.bedgraph



awk -v scaling_factor=2.31067063777344 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_iPSCpatient_IGG_R2.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_iPSCpatient_IGG_R2.dupmark.sorted.fragments.MG1655scaled.bedgraph



awk -v scaling_factor=4.18940754039497 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.fragments.MG1655scaled.bedgraph


awk -v scaling_factor=3.18445106295574 \
    -F'\t' -v OFS='\t' '{ $4 = $4 * scaling_factor; print $0 }' \
    output/bowtie2/8wN_iPSCpatient_IGG_R1.dupmark.sorted.fragments.bedgraph > output/bowtie2/8wN_iPSCpatient_IGG_R1.dupmark.sorted.fragments.MG1655scaled.bedgraph


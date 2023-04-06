#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=10G




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


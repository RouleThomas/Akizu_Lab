#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G
#SBATCH --time=50:00:00


x=("8wN_KO_IGG_R1" "8wN_KO_H3K27me3_R1"
   "8wN_KO_IGG_R2" "8wN_KO_H3K27me3_R2")




for x in "${x[@]}"; do
    bedtools bamtobed -bedpe -i output/bowtie2/${x}.dupmark.sorted.bam > output/bowtie2/${x}.dupmark.sorted.bed
    awk '$1==$4 && $6-$2 < 1000 {print $0}' output/bowtie2/${x}.dupmark.sorted.bed > output/bowtie2/${x}.dupmark.sorted.clean.bed # Filter out >1000bp fragment
    cut -f 1,2,6 output/bowtie2/${x}.dupmark.sorted.clean.bed | sort -k1,1 -k2,2n -k3,3n > output/bowtie2/${x}.dupmark.sorted.fragments.bed
    bedtools genomecov -bg -i output/bowtie2/${x}.dupmark.sorted.fragments.bed -g ../../Master/meta/GRCh38_chrom_sizes.tab > output/bowtie2/${x}.dupmark.sorted.fragments.bedgraph
done
  

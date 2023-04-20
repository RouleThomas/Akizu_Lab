#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=250G
#SBATCH --time=100:00:00


input_list=("2dN_HET_H3K27me3_R1"
"2dN_HET_H3K27me3_R2"
"2dN_HET_input_R1"
"2dN_HET_input_R2"
"2dN_KO_H3K27me3_R1"
"2dN_KO_H3K27me3_R2"
"2dN_KO_input_R1"
"2dN_KO_input_R2"
"2dN_WT_H3K27me3_R1"
"2dN_WT_H3K27me3_R2"
"2dN_WT_input_R1"
"2dN_WT_input_R2"
"NPC_HET_H3K27me3_R1"
"NPC_HET_H3K27me3_R2"
"NPC_HET_input_R1"
"NPC_HET_input_R2"
"NPC_KO_H3K27me3_R1"
"NPC_KO_H3K27me3_R2"
"NPC_KO_input_R1"
"NPC_KO_input_R2"
"NPC_WT_H3K27me3_R1"
"NPC_WT_H3K27me3_R2"
"NPC_WT_input_R1"
"NPC_WT_input_R2")

for x in "${input_list[@]}"; do
    bedtools bamtobed -i output/bowtie2_endtoend/${x}.dupmark.sorted.bam > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${x}.bed
done





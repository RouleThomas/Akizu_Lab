#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=250G
#SBATCH --time=100:00:00


input_list=("ESC_HET_H3K27me3_R1"
"ESC_HET_H3K27me3_R2"
"ESC_HET_input_R1"
"ESC_HET_input_R2"
"ESC_KO_H3K27me3_R1"
"ESC_KO_H3K27me3_R2"
"ESC_KO_input_R1"
"ESC_KO_input_R2"
"ESC_WT_H3K27me3_R1"
"ESC_WT_H3K27me3_R2"
"ESC_WT_H3K27me3_R3"
"ESC_WT_input_R1"
"ESC_WT_input_R2"
"ESC_WT_input_R3")

for x in "${input_list[@]}"; do
    bedtools bamtobed -i output/bowtie2_endtoend/${x}.dupmark.sorted.bam > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig/${x}.bed
done





#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G



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
    bamCoverage --bam output/bowtie2_endtoend/${x}.dupmark.sorted.bam \
        --outFileName output/bigwig/${x}.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 10 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done





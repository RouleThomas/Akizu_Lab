#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G



input_list=("NPC_HET_H3K27me3_R1"
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
    bamCoverage --bam output/bowtie2_endtoend/downsample/${x}.dupmark.sorted.bam \
        --outFileName output/bigwig_downsample/${x}.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 10 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done





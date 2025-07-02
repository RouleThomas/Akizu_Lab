#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00



input_list=(
"ESC_KO_EZH1"
"ESC_KO_EZH2"
"ESC_KO_H3K27me3"
"ESC_KO_H3K27me3epi"
"ESC_OEKO_EZH1"
"ESC_OEKO_EZH2"
"ESC_OEKO_H3K27me3"
"ESC_WT_EZH1"
"ESC_WT_EZH2"
"ESC_WT_H3K27me3"
"ESC_WT_H3K27me3epi"
"ESC_WT_IGG"
"NPC_KO_EZH1"
"NPC_KO_EZH2"
"NPC_KO_H3K27me3"
"NPC_KO_IGG"
"NPC_OEKO_EZH1"
"NPC_OEKO_EZH2"
"NPC_OEKO_H3K27me3"
"NPC_WT_EZH1"
"NPC_WT_EZH2"
"NPC_WT_H3K27me3"
)

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.unique.dupmark.sorted.bam \
        --outFileName output/bigwig/${x}.unique.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done


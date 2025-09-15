#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00



input_list=(
"ESC_KO_EZH1_R1"
"ESC_KO_EZH2_R1"
"ESC_KO_H3K27me3_R1"
"ESC_OEKO_EZH1_R1"
"ESC_OEKO_EZH2_R1"
"ESC_OEKO_H3K27me3_R1"
"ESC_WT_EZH1_R1"
"ESC_WT_EZH2_R1"
"ESC_WT_H3K27me3_R1"
"ESC_WT_IGG_R1"
"ESC_KO_EZH1_R2"
"ESC_KO_EZH2_R2"
"ESC_KO_H3K27me3_R2"
"ESC_OEKO_EZH1_R2"
"ESC_OEKO_EZH2_R2"
"ESC_OEKO_H3K27me3_R2"
"ESC_WT_EZH1_R2"
"ESC_WT_EZH2_R2"
"ESC_WT_H3K27me3_R2"
"ESC_WT_IGG_R2"
"ESC_KO_EZH1_R3"
"ESC_KO_EZH2_R3"
"ESC_KO_H3K27me3_R3"
"ESC_OEKO_EZH1_R3"
"ESC_OEKO_EZH2_R3"
"ESC_OEKO_H3K27me3_R3"
"ESC_WT_EZH1_R3"
"ESC_WT_EZH2_R3"
"ESC_WT_H3K27me3_R3"
"ESC_WT_IGG_R3"
)

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}_noXchr.unique.dupmark.sorted.bam \
        --outFileName output/bigwig/${x}_noXchr.unique.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done





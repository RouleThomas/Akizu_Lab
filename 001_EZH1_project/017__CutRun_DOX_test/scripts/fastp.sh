#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
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

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done





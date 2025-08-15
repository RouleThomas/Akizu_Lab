#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
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

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done





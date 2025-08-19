#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


macs2_out="output/macs2/broad_noIGG"
blacklist="../../Master/meta/hg38-blacklist.v2.bed"




samples=(
  # WT
  "ESC_WT_EZH1_R1" "ESC_WT_EZH1_R2" "ESC_WT_EZH1_R3"
  "ESC_WT_EZH2_R1" "ESC_WT_EZH2_R2" "ESC_WT_EZH2_R3"
  "ESC_WT_H3K27me3_R1" "ESC_WT_H3K27me3_R2" "ESC_WT_H3K27me3_R3"
  "ESC_WT_IGG_R1" "ESC_WT_IGG_R2" "ESC_WT_IGG_R3"

  # KO
  "ESC_KO_EZH1_R1" "ESC_KO_EZH1_R2" "ESC_KO_EZH1_R3"
  "ESC_KO_EZH2_R1" "ESC_KO_EZH2_R2" "ESC_KO_EZH2_R3"
  "ESC_KO_H3K27me3_R1" "ESC_KO_H3K27me3_R2" "ESC_KO_H3K27me3_R3"

  # OEKO
  "ESC_OEKO_EZH1_R1" "ESC_OEKO_EZH1_R2" "ESC_OEKO_EZH1_R3"
  "ESC_OEKO_EZH2_R1" "ESC_OEKO_EZH2_R2" "ESC_OEKO_EZH2_R3"
  "ESC_OEKO_H3K27me3_R1" "ESC_OEKO_H3K27me3_R2" "ESC_OEKO_H3K27me3_R3"
)

for sample in "${samples[@]}"; do
  macs2 callpeak -t "output/bowtie2/${sample}.unique.dupmark.sorted.bam" \
    -f BAMPE --keep-dup auto \
    --nomodel -g hs \
    --outdir "${macs2_out}" -n "${sample}" --broad
done





#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7


input_list=(
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





for x in "${input_list[@]}"; do
    bigWigToBedGraph output/bigwig/${x}_noXchr.unique.dupmark.sorted.bw output/bigwig/${x}_noXchr.unique.dupmark.sorted.bedGraph
done







